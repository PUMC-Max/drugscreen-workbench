"""单浓度、连续计数药筛。仅进行描述性分析，不推断跨实验显著性。"""
import math, re, statistics as st, json, sqlite3, hashlib
from pathlib import Path
from datetime import datetime
import qc
from plate_model import cell_name, mixed_cells, cell_views, analysis_records, pair_status

VERSION = '1.3.0'
ROLES = ['药物库', 'DMSO对照', '阳性药物对照', '仅细胞', '参考药物', '原库空孔', '未使用']
TIMES = list(range(0, 25, 2))

def ranked_treatment(well):
    """显式允许参考处理/阳性药物同时出榜，不改变其对照身份。"""
    return well['role']=='药物库' or (well['role'] in ('参考药物','阳性药物对照') and well.get('include_in_ranking',False))

def parse_txt(text):
    blocks = re.split(r'Time Stamp:', text)
    result = {}
    for block in blocks[1:]:
        m = re.search(r'Elapsed:\s*([\d.]+)\s*hours', block)
        if not m:
            raise ValueError('时间块缺少 Elapsed 小时数')
        t = float(m.group(1))
        if t in result:
            raise ValueError(f'重复时间点：{t:g} h')
        wells = {}
        for line in block.splitlines():
            if not re.match(r'^[A-H]\t', line):
                continue
            cells = line.split('\t')
            if len(cells) != 13:
                raise ValueError(f'{t:g} h 的 {cells[0]} 行不是 12 个数值')
            for col, value in enumerate(cells[1:], 1):
                key = f'{cells[0]}{col:02}'
                if key in wells:
                    raise ValueError(f'重复孔位：{t:g} h / {key}')
                try:
                    n = float(value.strip())
                except ValueError:
                    raise ValueError(f'无法读取计数：{t:g} h / {key}')
                if not math.isfinite(n) or n < 0:
                    raise ValueError(f'计数必须是有限非负数：{t:g} h / {key}')
                wells[key] = n
        if len(wells) != 96:
            raise ValueError(f'{t:g} h 只有 {len(wells)} 孔，应为 96 孔')
        result[t] = wells
    if sorted(result) != TIMES:
        raise ValueError('需要完整的 0、2、4……24 h 共 13 个时间点；缺失值不自动插补。')
    return {'header': blocks[0].strip(), 'times': {str(int(t)): result[t] for t in sorted(result)}}

def read_txt(path):
    raw = Path(path).read_bytes()
    for enc in ('utf-8-sig', 'utf-16', 'gb18030'):
        try:
            text = raw.decode(enc)
            if 'Time Stamp:' in text:
                return text, parse_txt(text)
        except UnicodeError:
            pass
    raise ValueError('无法识别文件编码或未找到 Time Stamp 时间块')

def mean(values):
    v = [x for x in values if x is not None and math.isfinite(x)]
    return st.mean(v) if v else None

def validate_layout(layout):
    expected = {f'{r}{c:02}' for r in 'ABCDEFGH' for c in range(1,13)}
    if len(layout) != 96 or {w['well'] for w in layout} != expected:
        raise ValueError('板位必须包含 96 个不重复的孔')
    for w in layout:
        if w['role'] not in ROLES:
            raise ValueError('未知孔类型')
        if w['role'] == '药物库' and (not w['drug_id'] or not w['name'] or w['dose'] != 5):
            raise ValueError(f"{w['well']}：药物库孔需有药物编号、名称，且浓度为 5 μM")
    if not any(w['role']=='DMSO对照' and not w.get('excluded') for w in layout):
        raise ValueError('至少需要一个未排除的 DMSO 对照孔')
    if mixed_cells(dict(layout=layout)):
        if any(cell_name(w) not in ('EV','R273H') for w in layout if w['role'] not in ('未使用','原库空孔')):
            raise ValueError('混合细胞板的含细胞孔需要明确细胞类型')
        for cell in ('EV','R273H'):
            if not any(cell_name(w)==cell and w['role']=='DMSO对照' and not w.get('excluded') for w in layout):
                raise ValueError(cell+' 至少需要一个未排除的 DMSO 对照孔')

def analyze(record, qc_settings=None):
    if mixed_cells(record):
        metrics={};drugs=[];warnings=[];folds={}
        for view in cell_views(record):
            cell=view['analysis_cell'];part=analyze(view,qc_settings)
            metrics.update({w['well']:part['metrics'][w['well']] for w in record['layout'] if cell_name(w)==cell})
            drugs.extend(dict(d,cell=cell) for d in part['drugs'])
            warnings.extend(cell+'：'+warning for warning in part['warnings']);folds[cell]=part['control_fold']
        return dict(metrics=metrics,drugs=drugs,warnings=warnings,control_fold=None,control_folds=folds,qc=qc.inspect(record,qc_settings))
    layout, times = record['layout'], record['data']['times']
    c = [w['well'] for w in layout if w['role']=='DMSO对照' and not w.get('excluded')]
    positive_c = [w for w in c if times['0'][w] > 0]
    warnings = []
    if len(c) < 3:
        warnings.append('有效 DMSO 对照少于 3 孔，Z 分数稳定性有限')
    if len(positive_c) != len(c):
        warnings.append('存在 0 h 计数为零的 DMSO 孔，其增长比不可计算')
    med = {t: st.median(times[t][w] for w in c) if c else 0 for t in times}
    fold = {t: st.median(times[t][w]/times['0'][w] for w in positive_c) if positive_c else 0 for t in times}
    if fold['24'] <= 1:
        warnings.append('DMSO 中位增长倍数≤1，GI、GR 不可计算')
    elif fold['24'] < 1.2:
        warnings.append('DMSO 增长较小（<1.2 倍），GI、GR 可能不稳定；这是提示阈值')
    control24 = [times['24'][w] for w in c]
    sd = st.stdev(control24) if len(c)>1 else 0
    mad = st.median(abs(x-med['24']) for x in control24) if c else 0
    if not sd or not mad:
        warnings.append('对照 SD 或 MAD 为零，相应 Z 分数不可计算')
    metrics, groups = {}, {}
    for w in layout:
        key = w['well']; n0 = times['0'][key]; n24 = times['24'][key]
        ratio = n24/n0 if n0>0 else None
        v = n24/med['24'] if med['24']>0 else None
        try:
            gr = (2**(math.log(ratio)/math.log(fold['24']))-1 if ratio else -1) if ratio is not None and fold['24']>1 else None
        except OverflowError:
            gr = None
            warnings.append(f'{key}：GR 数值溢出，请检查对照增长和计数')
        gi = 1-(ratio-1)/(fold['24']-1) if ratio is not None and fold['24']>1 else None
        curve = [times[str(t)][key]/med[str(t)] if med[str(t)]>0 else None for t in TIMES]
        auc = 1-sum((curve[i]+curve[i+1])/2*2 for i in range(12))/24 if all(x is not None for x in curve) else None
        z = (n24-st.mean(control24))/sd if sd else None
        rz = (n24-med['24'])/(1.4826*mad) if mad else None
        metrics[key] = dict(V=v, GI=gi, GR=gr, Z=z, RZ=rz, TI=auc, ratio=ratio)
        if ranked_treatment(w) and not w.get('excluded'):
            groups.setdefault(w['drug_id'], []).append(w)
    drugs = []
    for drug_id, wells in groups.items():
        m = {k: mean([metrics[w['well']][k] for w in wells]) for k in ('V','GI','GR','Z','RZ','TI')}
        # 聚合必须覆盖全部有效复孔，不用一个可算孔掩盖另一个失效孔。
        for k in m:
            if any(metrics[w['well']][k] is None for w in wells): m[k] = None
        flags=[]
        if len(wells)!=2: flags.append(f'有效复孔 {len(wells)} 个（预期 2）')
        if any(times['0'][w['well']]==0 for w in wells): flags.append('含 0 h 零计数孔')
        drugs.append(dict(id=drug_id,name=wells[0]['name'],wells=','.join(w['well'] for w in wells),n=len(wells),flags='；'.join(flags),**m))
    report=qc.inspect(record,qc_settings)
    for d in drugs:d['flags']='；'.join(filter(None,[d['flags'],*qc.summaries(report,d['wells'].split(','))]))
    return dict(metrics=metrics,drugs=drugs,warnings=warnings,control_fold=fold['24'],qc=report)

def rankings(records, ev_min=None, h_max=None, delta_min=None, weight=.5, qc_settings=None):
    by_plate={r['plate']:r for r in analysis_records(records)}
    rows=[]
    for plate,e in by_plate.items():
        if not plate.endswith('E') or plate[:-1]+'H' not in by_plate: continue
        h=by_plate[plate[:-1]+'H']
        ea,ha=analyze(e,qc_settings),analyze(h,qc_settings)
        hd={d['id']:d for d in ha['drugs']}
        for ev in ea['drugs']:
            if ev['id'] not in hd: continue
            hh=hd[ev['id']]
            dv=ev['V']-hh['V'] if ev['V'] is not None and hh['V'] is not None else None
            dg=ev['GR']-hh['GR'] if ev['GR'] is not None and hh['GR'] is not None else None
            if dv is None: continue
            if ev_min is not None and ev['V']<ev_min: continue
            if h_max is not None and hh['V']>h_max: continue
            if delta_min is not None and dv<delta_min: continue
            # 两个独立目标的几何效用：H 抑制强度和 EV 保留程度。仅为用户可调排序规则。
            clip=lambda v:max(0,min(1,v))
            a,b=clip(1-hh['V']),clip(ev['V'])
            score=100*(a**weight if weight else 1)*(b**(1-weight) if weight<1 else 1)
            rows.append(dict(pair=plate[:-1],id=ev['id'],name=ev['name'],EV=ev['V'],H=hh['V'],dV=dv,dGR=dg,
                EVGR=ev['GR'],EVGI=ev['GI'],EVZ=ev['Z'],EVRZ=ev['RZ'],EVTI=ev['TI'],
                HGI=hh['GI'],HGR=hh['GR'],HZ=hh['Z'],HRZ=hh['RZ'],HTI=hh['TI'],score=score,
                flags='；'.join(filter(None,[ev['flags'],hh['flags'],*ea['warnings'],*ha['warnings']])),))
    return rows

class Store:
    def __init__(self,path):
        self.path=Path(path); self.path.parent.mkdir(parents=True,exist_ok=True)
        self.db=sqlite3.connect(self.path)
        self.db.execute('CREATE TABLE IF NOT EXISTS records (id INTEGER PRIMARY KEY, plate TEXT, created TEXT, active INTEGER, payload TEXT)')
        self.db.execute('CREATE UNIQUE INDEX IF NOT EXISTS active_plate ON records(plate) WHERE active=1')
        self.db.commit()
    def records(self,active=True):
        query='SELECT id,plate,created,active,payload FROM records'+(' WHERE active=1' if active else '')+' ORDER BY id DESC'
        return [dict(json.loads(p),id=i,plate=b,created=c,active=a) for i,b,c,a,p in self.db.execute(query)]
    def backup(self):
        folder=self.path.parent/'backups'; folder.mkdir(exist_ok=True)
        target=sqlite3.connect(folder/(datetime.now().strftime('%Y%m%d-%H%M%S-%f')+'.sqlite'))
        self.db.backup(target); target.close()
    def save(self,record,replace=False):
        validate_layout(record['layout'])
        record=dict(record,qc_snapshot=qc.inspect(record,record.get('qc_settings')))
        duplicate=next((r for r in self.records() if r['plate']==record['plate']),None)
        if duplicate and not replace: raise ValueError('该板已有有效读板结果，请明确替换')
        signature=hashlib.sha256(record['raw'].encode()).hexdigest()
        if any(r.get('hash')==signature and r['plate']!=record['plate'] for r in self.records()):
            raise ValueError('同一份读板结果已用于其他板号，请检查，避免重复导入')
        self.backup()
        with self.db:
            if duplicate: self.db.execute('UPDATE records SET active=0 WHERE id=?',(duplicate['id'],))
            self.db.execute('INSERT INTO records(plate,created,active,payload) VALUES(?,?,1,?)',
                (record['plate'],datetime.now().isoformat(timespec='seconds'),json.dumps(dict(record,hash=signature,version=VERSION),ensure_ascii=False)))
    def remove(self,record_id):
        self.backup()
        with self.db:self.db.execute('UPDATE records SET active=0 WHERE id=?',(record_id,))
