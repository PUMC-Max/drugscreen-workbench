"""从 PlateMap 和独立 TXT 构建文件示例。Excel 的 Data 表从不参与计算。"""
import argparse
import copy
import hashlib
import html
import json
import re
from pathlib import Path
from collections import Counter
import openpyxl
from engine import TIMES, Store, read_txt, analyze, rankings, ranked_treatment, validate_layout
import qc
import research
import reviews
from demo import write_csv
from plate_model import cell_views

HERE = Path(__file__).resolve().parent
DEFAULT_MAP = HERE/'example'/'example-input.xlsx'
DEFAULT_COUNTS = HERE/'example'/'example-counts.txt'
DEFAULT_OUT = HERE/'example'/'file_example'


def map_rows(rows):
    """保留原孔位和处理身份；绘图开关不当作排除标志。"""
    rows=iter(rows);header=next(rows)
    required={'plate_id','well','cell_line','treatment','replicate','treatment_class','dmso_reference'}
    if not required.issubset(header):raise ValueError('PlateMap 缺少列：'+', '.join(sorted(required-set(header))))
    layout=[];plate_ids=set()
    for row in rows:
        d=dict(zip(header,row))
        if not any(v is not None for v in row):continue
        well=str(d.get('well') or '').strip().upper();match=re.fullmatch(r'([A-H])0?([1-9]|1[0-2])',well)
        if not match:raise ValueError('未知孔位：'+well)
        well=f'{match[1]}{int(match[2]):02}'
        cell=str(d.get('cell_line') or '').strip()
        if cell not in ('EV','R273H'):raise ValueError(well+'：需明确 EV 或 R273H')
        name=str(d.get('treatment') or '').strip();kind=d.get('treatment_class')
        if not name:raise ValueError(well+'：缺少处理名称')
        ref=d.get('dmso_reference')
        if ref not in (True,False,None):raise ValueError(well+'：dmso_reference 需要布尔值')
        if ref and name.upper()!='DMSO':raise ValueError(well+'：非 DMSO 处理标为 DMSO 参考，请核对')
        if name.upper()=='DMSO':
            if not ref:raise ValueError(well+'：DMSO 参考标志未启用，请核对')
            role='DMSO对照'
        elif name.lower()=='cells only':role='仅细胞'
        elif name.upper()=='RSL3':role='阳性药物对照'
        elif kind=='Screen drug':role='药物库'
        elif kind=='Named drug':role='参考药物'
        else:raise ValueError(well+'：无法确定处理类型 '+str(kind))
        dose=5 if role=='药物库' else 1 if role=='阳性药物对照' else None
        include=role in ('药物库','阳性药物对照','参考药物')
        layout.append(dict(well=well,cell=cell,role=role,name=name,drug_id=name if include else '',dose=dose,
            dose_source='既定实验设置，非本表字段' if dose is not None else '表中未提供',
            replicate=d.get('replicate'),excluded=False,reason='',include_in_ranking=include,
            source_row=d,source_include_in_plots=d.get('include_in_plots')))
        if d.get('plate_id'):plate_ids.add(str(d['plate_id']).strip())
    if len(plate_ids)>1:raise ValueError('该 TXT 示例仅对应一块物理板，PlateMap 含多个板号')
    validate_layout(layout)
    for cell in ('EV','R273H'):
        if not any(w['cell']==cell and w['role']=='DMSO对照' for w in layout):raise ValueError(cell+' 缺少自身 DMSO 对照')
    return layout,next(iter(plate_ids),None)


def load_source(map_path=DEFAULT_MAP,counts_path=DEFAULT_COUNTS):
    map_path=Path(map_path);counts_path=Path(counts_path)
    wb=openpyxl.load_workbook(map_path,read_only=True,data_only=True)
    try:
        if 'PlateMap' not in wb.sheetnames:raise ValueError('需要 PlateMap 工作表')
        layout,plate=map_rows(wb['PlateMap'].values)
    finally:wb.close()
    raw,data=read_txt(counts_path)
    return dict(plate=plate or counts_path.stem,plate_id_source='PlateMap' if plate else '板号为空，以 TXT 文件名作为示例标识',
        layout=layout,raw=raw,data=data,filename=counts_path.name,synthetic=False,example=True,
        provenance=dict(layout_file=str(map_path.resolve()),layout_sheet='PlateMap',layout_sha256=hashlib.sha256(map_path.read_bytes()).hexdigest(),
            counts_file=str(counts_path.resolve()),counts_sha256=hashlib.sha256(counts_path.read_bytes()).hexdigest(),
            note='只读取 PlateMap 板位；counts 全部来自 TXT；所有 22 个药物处理含 RSL3 入榜。'))


def analysis_views(source,source_id=1):
    """一块物理板的两个计算分组，保留原始坐标、原始计数及来源。"""
    return cell_views(dict(source,id=source_id,source_snapshot=source_id))


def save_source(store,source):
    # 物理板完整保存一次；E/H 仅是内存分析视图，不伪造两次独立实验。
    store.db.execute('CREATE TABLE IF NOT EXISTS file_example_sources (id INTEGER PRIMARY KEY, hash TEXT UNIQUE NOT NULL, created TEXT NOT NULL, payload TEXT NOT NULL)')
    payload=json.dumps(source,ensure_ascii=False,sort_keys=True);digest=hashlib.sha256(payload.encode()).hexdigest()
    found=store.db.execute('SELECT id FROM file_example_sources WHERE hash=?',(digest,)).fetchone()
    if found:return found[0]
    store.backup()
    with store.db:
        cur=store.db.execute('INSERT INTO file_example_sources(hash,created,payload) VALUES(?,?,?)',(digest,reviews.now(),payload))
    return cur.lastrowid


def natural(name):
    m=re.fullmatch(r'drug(\d+)',name)
    return (0,int(m[1])) if m else (1,name.lower())


def result_signature(records,result):
    return json.dumps(dict(records=[(r['id'],r['plate'],r['hash']) for r in records],baseline=result['baseline'],top_k=result['top_k'],qc=qc.config()),sort_keys=True)


def make_figures(out,source,views,profiles,stability):
    import plots,research_plots
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    fig=plots.Figure(figsize=(15,8),dpi=130);FigureCanvasAgg(fig);ax=fig.add_subplot();fig.subplots_adjust(left=.04,right=.99,top=.87,bottom=.12)
    colors={'DMSO对照':'#dce8ec','仅细胞':'#ebedef','阳性药物对照':'#f8dbb4'}
    for w in source['layout']:
        x=int(w['well'][1:])-1;y=ord(w['well'][0])-65
        fill=colors.get(w['role'],'#dceef7' if w['cell']=='EV' else '#eddef6')
        ax.add_patch(plots.Rectangle((x-.48,y-.45),.96,.9,facecolor=fill,edgecolor='white'))
        name='仅细胞' if w['role']=='仅细胞' else w['name']
        ax.text(x,y,f"{w['well']}\n{name}",ha='center',va='center',fontsize=9)
    ax.set(xlim=(-.6,11.6),ylim=(7.6,-.6),xticks=range(12),yticks=range(8));ax.set_xticklabels(range(1,13));ax.set_yticklabels(list('ABCDEFGH'));ax.xaxis.tick_top()
    ax.axvline(5.5,color='#7a8b94',lw=1.5);ax.set_title('PlateMap 原始板位 · 左侧 EV / 右侧 R273H',loc='left',fontsize=18,pad=28)
    fig.text(.04,.04,'一块物理板；每个处理在每种细胞中各两个技术复孔。板位来自 Excel，读数仅来自 TXT。',fontsize=11)
    fig.savefig(out/'plate_map.png')
    names=sorted({p['drug_id'] for p in profiles},key=natural)
    for page,start in enumerate(range(0,len(names),12),1):
        chosen=names[start:start+12];fig=plots.Figure(figsize=(14,10.4),dpi=130);FigureCanvasAgg(fig)
        fig.subplots_adjust(left=.06,right=.98,bottom=.10,top=.88,hspace=.72,wspace=.30)
        for i,drug in enumerate(chosen):
            ax=fig.add_subplot(3,4,i+1);vals=[]
            for cell,color in [('EV',plots.BLUE),('R273H',plots.PURPLE)]:
                p=next(p for p in profiles if p['drug_id']==drug and p['cell']==cell)
                for rep in p['replicates']:ax.plot(TIMES,rep['curve'],color=color,alpha=.4,ls='--',lw=1)
                ax.plot(TIMES,p['curve'],color=color,label=cell,lw=2);vals.extend(v for r in p['replicates'] for v in r['curve'] if v is not None)
            hp=next(p for p in profiles if p['drug_id']==drug and p['cell']=='R273H')
            ax.axhline(.8,color='#aa8156',ls=':',lw=1)
            ax.set(xlim=(0,24),ylim=(0,max([1.15]+vals)*1.05),xticks=[0,6,12,18,24],xlabel='Elapsed（h）',ylabel='相对增长 R(t)')
            ax.set_title(f'{drug}\nH：{hp["label"]}',loc='left',fontsize=10);ax.grid(alpha=.15)
            if i==0:ax.legend(fontsize=8)
        fig.suptitle(f'文件示例 · 时序反应分型（{page}）',x=.06,ha='left',fontsize=18)
        fig.text(.06,.035,'细线为技术复孔，粗线为均值；各自使用同细胞 DMSO 归一化。阈值仅供探索；保留异常及复孔差异。',fontsize=10)
        fig.savefig(out/f'temporal_{page}.png')
    fig=research_plots.stability_figure(stability,len(names));FigureCanvasAgg(fig);fig.savefig(out/'stability.png',dpi=130)


def generate(map_path=DEFAULT_MAP,counts_path=DEFAULT_COUNTS,destination=DEFAULT_OUT):
    out=Path(destination);out.mkdir(parents=True,exist_ok=True);workspace=out/'workspace';workspace.mkdir(exist_ok=True)
    source=load_source(map_path,counts_path);store=Store(workspace/'data'/'screen.sqlite')
    try:
        sid=save_source(store,source);views=analysis_views(source,sid)
        profiles=research.temporal_profiles(views);ranks=rankings(views);stab=research.stability(views,top_k=5)
        stab['source_signature']=result_signature(views,stab);stab['provenance']=source['provenance'];stab['qc_parameters']=qc.config()
        runs=reviews.runs(store,'stability')
        if not runs or runs[0]['result'].get('source_signature')!=stab['source_signature']:reviews.save_run(store,'stability',stab)
        events=reviews.history(store)
        reports={r['analysis_cell']:analyze(r) for r in views}
        findings=[dict(cell=cell,**f) for cell,a in reports.items() for f in a['qc']['findings']]
        payload=dict(example=True,synthetic=False,physical_plates=1,provenance=source['provenance'],source_snapshot=sid,
            source=source,analysis_views=[dict(id=r['id'],plate=r['plate'],physical_plate=r['physical_plate'],cell=r['analysis_cell'],hash=r['hash']) for r in views],
            temporal_parameters=research.temporal_config(),profiles=profiles,rankings=ranks,stability=stab,qc=findings,reviews=events)
        (out/'analysis.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
        for filename,rows in [('plate_map.csv',source['layout']),('rankings.csv',sorted(ranks,key=lambda r:-r['dV'])),('temporal_profiles.csv',profiles),
                              ('stability_summary.csv',stab['summary']),('stability_scenarios.csv',stab['details']),('qc.csv',findings),('rescreen_history.csv',events)]:write_csv(out/filename,rows)
        make_figures(out,source,views,profiles,stab)
        esc=lambda x:html.escape(str(x));fmt=lambda x:'不可计算' if x is None else f'{x:.3f}'
        hp={p['drug_id']:p for p in profiles if p['cell']=='R273H'};sp={s['drug_id']:s for s in stab['summary']}
        table=''.join('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in [r['name'],fmt(r['EV']),fmt(r['H']),fmt(r['dV']),fmt(r['dGR']),hp[r['id']]['label'],hp[r['id']]['agreement'],f"{sp[r['id']]['top_hits']}/{sp[r['id']]['total']}"])+'</tr>' for r in sorted(ranks,key=lambda r:-r['dV']))
        qc_table=''.join('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in [f['cell'],f['level'],f['title'],','.join(f['wells']),f['reason']])+'</tr>' for f in findings)
        controls='；'.join(f"{cell}：{', '.join(w['well'] for w in source['layout'] if w['cell']==cell and w['role']=='DMSO对照')}，24 h 中位增长 {a['control_fold']:.3f} 倍" for cell,a in reports.items())
        review_text='尚无人工复筛记录。可在“启动文件示例.cmd”的复筛记录页填写计划、证据及结果；系统不会自动生成已复现标签。' if not events else f'已有 {len(events)} 条人工复筛历史版本，详见记录导出或桌面示例。'
        pages=''.join(f'<img src="temporal_{i}.png" alt="时序分型第{i}页">' for i in range(1,3))
        report=f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>指定文件药筛示例</title>
<style>body{{font:16px/1.7 "Microsoft YaHei",sans-serif;color:#233943;background:#f1f5f6;margin:0}}main{{max-width:1250px;margin:28px auto;background:white;padding:32px;border-radius:16px}}h1,h2{{color:#14566b}}.notice{{padding:16px;background:#fff0d8;border-left:5px solid #d18a2b}}img{{width:100%;height:auto}}table{{width:100%;border-collapse:collapse;font-size:14px}}th,td{{text-align:left;padding:9px;border-bottom:1px solid #dbe3e8}}th{{background:#e9f1f4}}a{{color:#177ca0}}.scroll{{overflow-x:auto}}</style>
<main><h1>指定文件的药筛示例</h1><p class="notice">板位仅来自 <b>{esc(Path(map_path).name)} / PlateMap</b>；counts 仅来自 <b>{esc(Path(counts_path).name)}</b>。<br>一块物理板，EV / R273H 各48孔，22个药物处理、每种细胞各两个技术复孔、13个时间点。按用户要求作为示例，未并入正式库。全部结果均为这份TXT的实际计算值。</p>
<p>板位未填板号，以TXT文件名标识。drug1–drug16 沿用既定 5 μM 设置，RSL3 沿用既定 1 μM 设置；其余已命名处理浓度未提供，保持未填写。全部药物处理含RSL3参与排名，DMSO和仅细胞不参与药物排名。RSL3仍保留阳性对照身份。</p>
<h2>板位与归一化参考</h2><img src="plate_map.png" alt="96孔板位"><p>{esc(controls)}。每类细胞只使用自身DMSO参考；两个技术复孔不算独立实验。</p>
<h2>全部22项结果</h2><p>按ΔV由高到低展示全部处理，保留负值和质控提示。EV/H列为24 h的DMSO归一化V。稳定性Top-5沿用A类候选定义：ΔV&gt;0 且 H V&lt;1；未满足条件会显示筛除次数，原始指标仍保留。</p><div class="scroll"><table><tr><th>处理</th><th>EV V</th><th>H V</th><th>ΔV</th><th>ΔGR</th><th>H时序类别</th><th>H复孔分类</th><th>Top-5次数</th></tr>{table}</table></div>
<h2>时序反应分型</h2><p>R(t)=(本孔N(t)/N(0))/同细胞DMSO增长中位数。连续3点≤0.8触发持续抑制判据，迟发起点为12 h；这些是可调整的探索性参数。</p>{pages}
<h2>排名稳定性</h2><p>{len(stab['scenarios'])}种分析设置，比较算法、门槛、权重及单复孔选择。Top-5频率的分母为全部设置，不能解释为复筛成功概率。</p><img src="stability.png" alt="排名稳定性">
<h2>复筛记录</h2><p>{review_text}</p><h2>异常与检查提示</h2><p>下列提示不删改数据、不自动排除，不会因阈值触发而将可计算结果置为NA。</p><div class="scroll"><table><tr><th>细胞</th><th>级别</th><th>检查项</th><th>孔位</th><th>说明</th></tr>{qc_table}</table></div>
<h2>数据导出</h2><p><a href="plate_map.csv">板位</a> · <a href="rankings.csv">全部药效指标</a> · <a href="temporal_profiles.csv">时序分型</a> · <a href="stability_summary.csv">稳定性汇总</a> · <a href="stability_scenarios.csv">逐方案结果</a> · <a href="qc.csv">检查清单</a> · <a href="analysis.json">完整分析与来源</a></p>
<p>示例源数据快照 #{sid}。原始文件未修改；按仪器Elapsed计算0–24 h。复筛记录由人工填写，尚未训练预测模型。</p></main></html>'''
        (out/'示例报告.html').write_text(report,encoding='utf-8')
        return out,source,views,payload
    finally:store.db.close()


def prepare_workspace(map_path=DEFAULT_MAP,counts_path=DEFAULT_COUNTS,destination=DEFAULT_OUT):
    """只在第一次启动时导入；再次打开不覆盖修改，也不恢复已移除记录。"""
    workspace=Path(destination)/'workspace';workspace.mkdir(parents=True,exist_ok=True)
    store=Store(workspace/'data'/'screen.sqlite')
    try:
        history=store.records(False)
        if not history:
            source=load_source(map_path,counts_path)
            source['note']='指定文件示例：PlateMap板位 + TXT读数；22项处理全部参与分析'
            store.save(source);history=store.records(False)
        templates=workspace/'templates.json'
        if not templates.exists():
            templates.write_text(json.dumps({r['plate']:r['layout'] for r in reversed(history)},ensure_ascii=False),encoding='utf-8')
        return workspace
    finally:store.db.close()


def open_workspace(root,workspace):
    from app import App
    app=App(root,demo=True,workspace=workspace,example_label='指定文件示例')
    if app.research.result:app.research.top.set(str(app.research.result['top_k']));app.research.mark_stability()
    app.nb.select(app.vistab)
    return app


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--open',action='store_true');parser.add_argument('--smoke',action='store_true')
    parser.add_argument('--map',default=str(DEFAULT_MAP));parser.add_argument('--counts',default=str(DEFAULT_COUNTS));parser.add_argument('--out',default=str(DEFAULT_OUT));args=parser.parse_args()
    if args.open:
        import tkinter as tk
        workspace=prepare_workspace(args.map,args.counts,args.out)
        root=tk.Tk();app=open_workspace(root,workspace)
        if args.smoke:root.after(700,app.close)
        root.mainloop()
    else:
        out,source,views,payload=generate(args.map,args.counts,args.out);print(out)
