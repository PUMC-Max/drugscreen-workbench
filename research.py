"""AI4S 前置分析：可解释规则、确定性敏感性分析、人工复筛标签。"""
import copy, math, statistics as st
from engine import TIMES, rankings, VERSION, ranked_treatment
from plate_model import analysis_records
import qc

TEMPORAL_DEFAULTS=dict(inhibition=.8,recovery=.9,rebound=.25,late_hour=12,consecutive=3)
TEMPORAL_HELP='''时序分型是可解释的探索性规则，不是训练好的AI或机制结论。
R(t)=(本孔N(t)/N(0))/(有效DMSO各孔N(t)/N(0)的中位数)，再取有效复孔均值。
默认连续3个观测点R≤0.8，定义出现持续抑制；起点是这一段的第一个采样时间，时间分辨率2 h。
末3点均值≥0.9，且比曾满足持续抑制的最低点回升≥0.25：抑制后恢复。
起点≥12 h且末3点均值≤0.8：迟发抑制；早于12 h、末端仍受抑制、起点后至少75%点低于阈值：持续抑制。
未满足连续点要求但曾低于阈值：短暂/波动；始终未达阈值：未见持续抑制；其余：复杂/混合。
同时列出每个复孔的分类，一致不等于跨实验可重复。异常只提示，不屏蔽数值、不改变类别。
0 h为上机首拍，不是精确加药时刻。阈值可调，不能推断药物靶点、死亡或真实耐药。
'''
STABILITY_HELP='''排名稳定性是预设分析条件的敏感性比较，不是概率、P值或置信区间。
包括当前基准、ΔV/ΔGR/时间平均抑制差、综合权重0.25/0.50/0.75、门槛扰动与4种单复孔组合。
全部方案逐项展示；Top-K入选频率的分母为全部方案，不因筛除或指标不可计算而缩小。
同时报告可排序次数、筛除次数和不可计算次数。名次范围仅基于可排序方案。
并列采用竞赛排名（1,1,3），Top-K包含边界并列项，因此实际入选数可能超过K。
单复孔方案是对两个技术复孔的敏感性检查，不把它们当作独立实验，也不修改原始记录。
门槛比较沿用当前A/B榜定义。未设置门槛时用展示出的示例网格；不是经过验证的筛选标准。
'''

def temporal_config(values=None):
    p=dict(TEMPORAL_DEFAULTS);p.update(values or {})
    for k in TEMPORAL_DEFAULTS:
        p[k]=float(p[k])
        if not math.isfinite(p[k]):raise ValueError('时序阈值需要有限数字')
    if not 0<p['inhibition']<p['recovery'] or p['rebound']<=0:raise ValueError('要求 0<抑制阈值<恢复阈值，恢复幅度>0')
    if not 0<=p['late_hour']<=24 or p['consecutive']!=int(p['consecutive']) or not 2<=p['consecutive']<=6:raise ValueError('迟发起点需在0–24 h，连续点数需为2–6的整数')
    p['consecutive']=int(p['consecutive']);return p

def classify(values,settings=None):
    p=temporal_config(settings)
    if len(values)!=13 or any(v is None or not math.isfinite(v) for v in values):
        return dict(label='无法分型',onset=None,end=None,minimum=None,rebound=None,TI=None,reason='缺少可计算的增长比；未对原始药效结果做屏蔽')
    n=p['consecutive'];starts=[i for i in range(1,14-n) if all(v<=p['inhibition'] for v in values[i:i+n])]
    start=starts[0] if starts else None;end=st.mean(values[-3:]);minimum=min(values[1:]);reb=end-minimum
    if start is None:label='短暂/波动' if minimum<=p['inhibition'] else '未见持续抑制'
    elif end>=p['recovery'] and reb>=p['rebound']:label='抑制后恢复'
    elif TIMES[start]>=p['late_hour'] and end<=p['inhibition']:label='迟发抑制'
    elif end<=p['inhibition'] and sum(v<=p['inhibition'] for v in values[start:])/len(values[start:])>=.75:label='持续抑制'
    else:label='复杂/混合'
    ti=1-sum((values[i]+values[i+1]) for i in range(12))/24
    return dict(label=label,onset=TIMES[start] if start is not None else None,end=end,minimum=minimum,rebound=reb,TI=ti,
        reason=f'连续{n}点阈值≤{p["inhibition"]:g}；末3点均值={end:.3f}，最低={minimum:.3f}，回升={reb:.3f}')

def temporal_profiles(records,settings=None,qc_settings=None):
    p=temporal_config(settings);result=[]
    for record in analysis_records(records):
        ts=record['data']['times'];active=[w for w in record['layout'] if not w.get('excluded')]
        controls=[w['well'] for w in active if w['role']=='DMSO对照' and ts['0'][w['well']]>0]
        cf=[st.median(ts[str(t)][k]/ts['0'][k] for k in controls) if controls else 0 for t in TIMES]
        groups={}
        for w in active:
            if ranked_treatment(w):groups.setdefault(w['drug_id'],[]).append(w)
        report=qc.inspect(record,qc_settings)
        for drug,wells in groups.items():
            trajectories=[]
            for w in wells:
                k=w['well'];n0=ts['0'][k]
                curve=[ts[str(t)][k]/n0/c if n0>0 and c>0 else None for t,c in zip(TIMES,cf)]
                trajectories.append(dict(well=k,curve=curve,**classify(curve,p)))
            avg=[st.mean(r['curve'][i] for r in trajectories) if all(r['curve'][i] is not None for r in trajectories) else None for i in range(13)]
            labels=[r['label'] for r in trajectories];agreement='一致' if len(set(labels))==1 and len(labels)==2 else '仅一个有效复孔' if len(labels)==1 else '不一致 / 数量非预期'
            result.append(dict(plate=record['plate'],pair=record['plate'][:-1],record_id=record.get('id'),drug_id=drug,name=wells[0]['name'],
                cell='EV' if record['plate'].endswith('E') else 'R273H',curve=avg,replicates=trajectories,agreement=agreement,
                flags='；'.join(qc.summaries(report,[w['well'] for w in wells])),**classify(avg,p)))
    return result

def replicate_subset(records,eidx,hidx):
    out=copy.deepcopy(analysis_records(records))
    for r in out:
        groups={}
        for w in r['layout']:
            if ranked_treatment(w) and not w.get('excluded'):groups.setdefault(w['drug_id'],[]).append(w)
        for wells in groups.values():
            if len(wells)==2:
                selected=sorted(wells,key=lambda w:w['well'])[eidx if r['plate'].endswith('E') else hidx]['well']
                for w in wells:
                    if w['well']!=selected:w.update(excluded=True,reason='仅用于敏感性方案；原记录未改动')
    return out

def score(row,metric,weight):
    if metric=='ΔV':return row['dV']
    if metric=='ΔGR':return row['dGR']
    if metric=='时间抑制差':return row['HTI']-row['EVTI'] if row['HTI'] is not None and row['EVTI'] is not None else None
    if metric=='综合分':
        a=max(0,min(1,1-row['H']));b=max(0,min(1,row['EV']))
        return 100*(a**weight if weight else 1)*(b**(1-weight) if weight<1 else 1)
    mapping={'H GI':('HGI',1),'H GR':('HGR',-1),'H Z':('HZ',-1),'H Robust Z':('HRZ',-1),'H 时间平均抑制':('HTI',1)}
    if metric not in mapping:raise ValueError('不支持的排序指标')
    key,sign=mapping[metric];return row[key]*sign if row[key] is not None else None

def stability(records,baseline=None,top_k=10,qc_settings=None):
    base=dict(mode='A',metric='ΔV',weight=.5,ev_min=None,h_max=None,delta_min=None);base.update(baseline or {})
    if not isinstance(top_k,int) or not 1<=top_k<=100:raise ValueError('Top K应为1–100整数')
    base['weight']=float(base['weight'])
    if not 0<=base['weight']<=1:raise ValueError('权重需在0–1内')
    for key in ('ev_min','h_max','delta_min'):
        if base[key] in ('',None):base[key]=None
        else:
            base[key]=float(base[key])
            if not math.isfinite(base[key]):raise ValueError('门槛必须为有限数字')
    scenarios=[dict(base,name='当前基准',replicate=None)]
    for metric in ('ΔV','ΔGR','时间抑制差'):scenarios.append(dict(base,name='指标：'+metric,metric=metric,replicate=None))
    for w in (.25,.5,.75):scenarios.append(dict(base,name=f'综合权重{w:g}',metric='综合分',weight=w,replicate=None))
    for key,label,defaults in [('ev_min','EV门槛',(.7,.8,.9)),('h_max','H门槛',(.5,.7,.9)),('delta_min','ΔV门槛',(.1,.2,.3))]:
        choices=defaults if base[key] is None else tuple(base[key]+delta for delta in (-.1,0,.1))
        for v in choices:scenarios.append(dict(base,**{key:v},name=f'{label}{v:.2f}',replicate=None))
    for e,h in ((0,0),(0,1),(1,0),(1,1)):scenarios.append(dict(base,name=f'单复孔 E{e+1}/H{h+1}',replicate=(e,h)))
    # 去除完全相同的方案，避免重复基准人为放大频率。
    unique=[];seen=set()
    for s in scenarios:
        signature=tuple((k,str(s[k])) for k in sorted(s) if k!='name')
        if signature not in seen:unique.append(s);seen.add(signature)
    cached={None:rankings(records,qc_settings=qc_settings)}
    for e,h in ((0,0),(0,1),(1,0),(1,1)):cached[(e,h)]=rankings(replicate_subset(records,e,h),qc_settings=qc_settings)
    universe={(r['pair'],r['id']):r for r in cached[None]};details=[];by_drug={key:[] for key in universe}
    for si,s in enumerate(unique):
        rows={(r['pair'],r['id']):r for r in cached[s['replicate']]};eligible={};states={}
        for key in universe:
            row=rows.get(key);value=score(row,s['metric'],s['weight']) if row else None
            if value is None or not math.isfinite(value):states[key]='不可计算';continue
            ok=row['dV']>0 and (not str(s['mode']).startswith('A') or row['H']<1)
            ok=ok and all(s[k] is None or (row[rkey]>=s[k] if k!='h_max' else row[rkey]<=s[k]) for k,rkey in [('ev_min','EV'),('h_max','H'),('delta_min','dV')])
            if not ok:states[key]='筛除';continue
            eligible[key]=value;states[key]='可排序'
        # 排序一次生成竞赛名次，避免全库筛选时逐药物进行平方复杂度比较。
        rankmap={};groupvalue=None;grouprank=0
        for position,(key,value) in enumerate(sorted(eligible.items(),key=lambda item:-item[1]),1):
            if groupvalue is None or not math.isclose(value,groupvalue,rel_tol=1e-9,abs_tol=1e-12):groupvalue=value;grouprank=position
            rankmap[key]=grouprank
        for key,original in universe.items():
            rank=rankmap.get(key)
            item=dict(pair=key[0],drug_id=key[1],name=original['name'],scenario=si,scenario_name=s['name'],state=states[key],rank=rank,top=rank is not None and rank<=top_k)
            details.append(item);by_drug[key].append(item)
    summary=[]
    for key,r in universe.items():
        ds=by_drug[key];ranks=[d['rank'] for d in ds if d['rank'] is not None];hits=sum(d['top'] for d in ds)
        summary.append(dict(pair=key[0],drug_id=key[1],name=r['name'],top_hits=hits,total=len(unique),frequency=hits/len(unique),
            rank_min=min(ranks) if ranks else None,rank_median=st.median(ranks) if ranks else None,rank_max=max(ranks) if ranks else None,
            eligible=len(ranks),filtered=sum(d['state']=='筛除' for d in ds),unavailable=sum(d['state']=='不可计算' for d in ds),flags=r['flags']))
    summary.sort(key=lambda r:(-r['frequency'],r['rank_median'] if r['rank_median'] is not None else float('inf'),r['pair'],r['drug_id']))
    serial_scenarios=[dict(s,replicate=list(s['replicate']) if s['replicate'] is not None else None) for s in unique]
    return dict(version=VERSION,top_k=top_k,baseline=base,scenarios=serial_scenarios,summary=summary,details=details,
        records=[dict(id=r.get('id'),plate=r['plate'],hash=r.get('hash')) for r in records],note=STABILITY_HELP)
