"""可解释的探索性质控规则；只标记，不修改任何计数或自动排除孔。"""
import math
import statistics as st
from plate_model import mixed_cells,cell_views

QC_VERSION = '1.0'
TIMES = list(range(0, 25, 2))
DEFAULTS = dict(baseline_rel=.35, robust_z=3.5, spike_rel=.4, neighbor_rel=.25,
                jump_rel=.75, pair_rms=.30, pair_cv=.30, control_cv=.25,
                control_rel=.35, min_growth=1.2, positive_max=.8, spatial_rel=.30)
LABELS = {
    'baseline_rel':'起始计数最小相对偏差', 'robust_z':'稳健偏离倍数',
    'spike_rel':'局部尖峰/凹陷相对偏差', 'neighbor_rel':'尖峰两侧允许相对差',
    'jump_rel':'相邻时间点变化比例', 'pair_rms':'复孔曲线相对 RMS 差异',
    'pair_cv':'24 h 复孔 CV', 'control_cv':'DMSO 孔间 CV',
    'control_rel':'DMSO 单孔最小相对偏差', 'min_growth':'DMSO 最低增长倍数',
    'positive_max':'阳性对照最高相对增长（可留空）', 'spatial_rel':'起始行/列/边缘中位相对偏差'}
HELP = '''异常检查的阈值是可调整的探索性提示，不是已验证的验收标准。
各项比例用小数表示，例如 0.30 = 30%。没有提示不等于整板合格。
起始计数：与有效含细胞孔的 0 h 中位数比较，需同时超过相对偏差和稳健倍数；MAD=0 时只用相对偏差。
尖峰/凹陷：中间点偏离左右点均值，且左右点相近；不插值替换原值。
相邻跳变：绝对变化/ max(前一点计数, 本孔0 h计数×0.1, 1)。可能是真实快速药效。
复孔曲线：各孔除以自身0 h，计算两曲线 RMS 差异/平均曲线 RMS，另查24 h CV；不判断哪一个孔有错。
DMSO：每个时间点检查原始计数CV、单孔相对其余DMSO中位数的偏离；排除孔不作参考。
阳性对照：按名称和浓度分组，比较“本孔N24/N0 / DMSO中位N24/N0”，高于设置上限仅提示反应不足。
阳性上限留空时仅展示数值，不判定反应是否达标。该判定与药物剂量无关，不混合不同阳性处理。
空间分布：只检查药物库孔的0 h，按行、列、边缘与其余药物孔比较；各侧至少4孔，至少70%组内孔同向偏移。
此布局不同孔处理不同，空间提示不能证明边缘效应；不使用24 h药物效应自动校正空间偏差。
仅细胞孔和阳性孔不会加入DMSO参考。空孔/未使用孔不参加生物学异常检查。
已经手工排除的孔保留在清单中并显示排除原因，不参与质控参考或药效聚合。
异常仅提示：触发阈值不会将任何可计算结果变为NA，不自动删改原值或影响排名。只有人工排除才改变药效。
质控“无法评估”仅指某项检查缺少参考，不会因此将药效结果变为NA。数学无定义或计算模型不适用仍单独标注不可计算。
当前设置重查历史不会覆盖保存时的质控快照。
资料：https://www.ncbi.nlm.nih.gov/books/NBK83783/
'''

def config(values=None):
    out=dict(DEFAULTS)
    for k,v in (values or {}).items():
        if k not in DEFAULTS:continue
        if k=='positive_max' and (v is None or v==''):out[k]=None;continue
        try:v=float(v)
        except (TypeError,ValueError):raise ValueError(LABELS[k]+'需要数字')
        if not math.isfinite(v) or v<=0:raise ValueError(LABELS[k]+'必须为有限正数')
        out[k]=v
    return out

def deviation(value, references, relative, z):
    if len(references)<3:return False
    median=st.median(references)
    scale=1.4826*st.median(abs(v-median) for v in references)
    return abs(value-median)>relative*max(abs(median),1) and (scale==0 or abs(value-median)>z*scale)

def inspect(record, settings=None):
    if mixed_cells(record):
        findings=[]
        for view in cell_views(record):
            cell=view['analysis_cell']
            findings.extend(dict(f,cell=cell,reason=cell+'：'+f['reason']) for f in inspect(view,settings)['findings'])
        return dict(version=QC_VERSION,parameters=config(settings),findings=findings,
                    alerts=sum(f['level']=='需核对' for f in findings),unavailable=sum(f['level']=='无法评估' for f in findings))
    p=config(settings); ts=record['data']['times']; findings=[]
    active=[w for w in record['layout'] if w['role'] not in ('未使用','原库空孔') and not w.get('excluded')]
    def add(code,title,wells,times,reason,level='需核对'):
        findings.append(dict(code=code,title=title,wells=sorted(wells),times=list(times),reason=reason,level=level))
    def values(w):return [ts[str(t)][w] for t in TIMES]
    for w in record['layout']:
        if w.get('excluded'):add('excluded','人工排除',[w['well']],[],w.get('reason') or '未记录原因','已排除')
    baseline=[ts['0'][w['well']] for w in active]
    if len(baseline)<4:add('coverage','起始计数检查不足',[],[0],'有效含细胞孔不足4个，无法稳健比较起始分布','无法评估')
    for w in active:
        key=w['well'];v=values(key)
        if v[0]==0:add('zero_start','起始计数为零',[key],[0],'GI/GR与本孔增长曲线比较不可计算')
        elif deviation(v[0],baseline,p['baseline_rel'],p['robust_z']):
            add('baseline','起始计数偏离',[key],[0],f'本孔 {v[0]:.2f}；含细胞孔中位数 {st.median(baseline):.2f}。可能为铺板或早期药效差异')
        spikes=[];jumps=[]
        for i in range(1,12):
            neighbors=(v[i-1]+v[i+1])/2;den=max(neighbors,1)
            if abs(v[i]-neighbors)/den>p['spike_rel'] and abs(v[i-1]-v[i+1])/den<=p['neighbor_rel']:spikes.append(TIMES[i])
        for i in range(1,13):
            if abs(v[i]-v[i-1])/max(v[i-1],v[0]*.1,1)>p['jump_rel']:jumps.append(TIMES[i])
        if spikes:add('spike','孤立尖峰 / 凹陷',[key],spikes,f'偏离相邻均值 > {p["spike_rel"]:.0%}，且两侧相近；请核查对应图像')
        if jumps:add('jump','相邻时间点跳变',[key],jumps,f'相对前一点变化 > {p["jump_rel"]:.0%}；时间标记为区间终点，可能是真实药效')
    groups={}
    for w in active:
        if w['role']=='药物库' or (w['role'] in ('参考药物','阳性药物对照') and w.get('include_in_ranking',False)):
            groups.setdefault(w['drug_id'],[]).append(w['well'])
    for drug,keys in groups.items():
        if len(keys)!=2:
            add('replicate_count','有效复孔数异常',keys,[],f'{drug}：有效复孔 {len(keys)}，预期2');continue
        end=[ts['24'][w] for w in keys]
        if st.mean(end)>0:
            cv=st.stdev(end)/st.mean(end)
            if cv>p['pair_cv']:add('replicate_end','复孔终点不一致',keys,[24],f'{drug}：24 h CV={cv:.1%}，阈值 {p["pair_cv"]:.0%}')
        if any(ts['0'][w]<=0 for w in keys):
            add('replicate_unavailable','复孔曲线无法比较',keys,[0],'某复孔起始计数为零','无法评估');continue
        a,b=[[x/ts['0'][w] for x in values(w)] for w in keys]
        diff=math.sqrt(st.mean((x-y)**2 for x,y in zip(a,b)))
        denom=math.sqrt(st.mean(((x+y)/2)**2 for x,y in zip(a,b)))
        if denom>0 and diff/denom>p['pair_rms']:
            affected=[t for t,x,y in zip(TIMES,a,b) if abs(x-y)>p['pair_rms']*max((x+y)/2,.1)]
            add('replicate_curve','复孔全程曲线不一致',keys,affected,f'{drug}：归一起始值后相对 RMS 差异={diff/denom:.1%}。无法仅凭两孔确定错误孔')
    controls=[w['well'] for w in active if w['role']=='DMSO对照']
    if len(controls)<3:add('controls_few','DMSO 对照不足',controls,[],f'有效对照 {len(controls)} 孔，孔间离群检查至少需4孔','无法评估')
    spread=[];cvmax=0
    for t in TIMES:
        cvvalues=[ts[str(t)][k] for k in controls]
        if len(cvvalues)>=2 and st.mean(cvvalues)>0:
            cv=st.stdev(cvvalues)/st.mean(cvvalues);cvmax=max(cvmax,cv)
            if cv>p['control_cv']:spread.append(t)
    if spread:add('control_cv','DMSO 孔间波动较大',controls,spread,f'最大CV={cvmax:.1%}，阈值 {p["control_cv"]:.0%}；将影响归一化')
    if len(controls)==3:add('outlier_coverage','DMSO 离群检查不足',controls,[],'留一孔检查至少需要其余3个对照（总计4孔）','无法评估')
    for key in controls:
        bad=[t for t in TIMES if deviation(ts[str(t)][key],[ts[str(t)][k] for k in controls if k!=key],p['control_rel'],p['robust_z'])]
        if bad:add('control_outlier','DMSO 单孔偏离', [key],bad,'与其余有效DMSO孔比较，超过相对偏差和稳健阈值；请先核查再决定是否排除')
    folds=[ts['24'][k]/ts['0'][k] for k in controls if ts['0'][k]>0]
    cf=st.median(folds) if folds else None
    if cf is None:add('control_growth_unknown','DMSO 增长无法评估',controls,[0,24],'无有效起始值','无法评估')
    elif cf<p['min_growth']:add('control_growth','DMSO 增长不足',controls,[0,24],f'增长倍数 {cf:.3f} < {p["min_growth"]:g}；GI/GR可能不稳定')
    positives={}
    for w in active:
        if w['role']=='阳性药物对照':positives.setdefault((w['name'],w['dose']),[]).append(w['well'])
    if not positives:add('positive_missing','无有效阳性对照',[],[],'无法检查阳性反应','无法评估')
    for (name,dose),keys in positives.items():
        pf=[ts['24'][k]/ts['0'][k] for k in keys if ts['0'][k]>0]
        if not pf or cf is None or cf<=0:
            add('positive_unknown','阳性反应无法评估',keys,[0,24],f'{name} {dose} μM：无可用增长比','无法评估');continue
        ratio=st.median(pf)/cf
        if len(pf)<len(keys):add('positive_partial','部分阳性孔无有效起始值',keys,[0],'该组反应只覆盖有效起始值的阳性孔','无法评估')
        reason=f'{name} {dose} μM：相对增长={ratio:.3f}（阳性增长倍数/DMSO增长倍数）'
        if p['positive_max'] is not None and ratio>p['positive_max']:
            add('positive_weak','阳性反应未达设置要求',keys,[0,24],reason+f' > {p["positive_max"]:g}；需结合细胞对RSL3的敏感性判断')
        else:add('positive_observed','阳性反应观测值',keys,[0,24],reason+'；不代表整板通过验收','信息')
    drugs=[w['well'] for w in active if w['role']=='药物库']
    tested=0
    def spatial(label,keys):
        nonlocal tested
        rest=[k for k in drugs if k not in keys]
        if len(keys)<4 or len(rest)<4:return
        tested+=1
        m=st.median(ts['0'][k] for k in rest);v=st.median(ts['0'][k] for k in keys)
        rel=(v-m)/max(m,1)
        same=sum((ts['0'][k]-m)*(1 if rel>0 else -1)>p['spatial_rel']*max(m,1) for k in keys)/len(keys)
        if abs(rel)>p['spatial_rel'] and same>=.7:
            add('spatial','起始空间分布偏差',keys,[0],f'{label}：相对其余药物孔中位偏差 {rel:+.1%}，同向孔占比 {same:.0%}；提示不能证明边缘效应')
    for row in 'ABCDEFGH':spatial(row+'行',[k for k in drugs if k[0]==row])
    for col in range(1,13):spatial(str(col)+'列',[k for k in drugs if int(k[1:])==col])
    spatial('边缘',[k for k in drugs if k[0] in 'AH' or int(k[1:]) in (1,12)])
    if not tested:add('spatial_unknown','空间分布无法评估',[],[0],'分组或参考药物孔不足（各需至少4孔）','无法评估')
    return dict(version=QC_VERSION,parameters=p,findings=findings,
                alerts=sum(f['level']=='需核对' for f in findings),
                unavailable=sum(f['level']=='无法评估' for f in findings))

def summaries(report, wells=None):
    """药物关联本孔和整板对照提示；不把其他药物异常传播给所有药物。"""
    controls={'control_cv','controls_few','control_outlier','control_growth','control_growth_unknown','positive_missing','positive_unknown','positive_weak'}
    return sorted({f['title'] for f in report['findings'] if f['level'] in ('需核对','无法评估') and
                   (wells is None or f['code'] in controls or set(wells)&set(f['wells']))})
