"""确定性的公开友好合成示例，不读取任何用户实验文件。"""
from pathlib import Path
import json,csv,html,math
from engine import Store,TIMES
import research,reviews

NAMES=['持续选择性抑制','迟发选择性抑制','短暂抑制后恢复','双侧强抑制','复孔分歧','基本无反应','单点凹陷','温和迟发抑制']
CURVES=[
 [1,.9,.75,.6,.5,.45,.4,.35,.32,.30,.28,.26,.25],
 [1,1,.99,.99,.97,.93,.75,.55,.42,.35,.31,.29,.28],
 [1,.85,.65,.5,.45,.5,.6,.75,.85,.95,1,1.02,1.02],
 [1,.8,.65,.5,.4,.3,.25,.2,.18,.16,.14,.13,.12],
 [1,.93,.85,.78,.7,.63,.58,.55,.53,.52,.51,.5,.5],
 [1,1,1,.99,.99,.99,.98,.98,.98,.98,.98,.97,.97],
 [1,1,1,.99,.99,.99,.30,.99,.99,.99,.99,.99,.99],
 [1,.98,.94,.9,.87,.84,.82,.79,.77,.75,.73,.72,.7]]

def synthetic_records():
    records=[];positions=[(row,c) for row in 'AB' for c in range(2,12,2)][:8]
    locations={f'{row}{c+j:02}':(i,j) for i,(row,c) in enumerate(positions) for j in (0,1)}
    for cell in ('E','H'):
        layout=[];times={str(t):{} for t in TIMES}
        for row in 'ABCDEFGH':
            for col in range(1,13):
                key=f'{row}{col:02}';role='未使用';drug=name='';dose=0
                if col in (1,12):role='DMSO对照' if row in 'ABC' else '阳性药物对照' if row in 'DEF' else '仅细胞'
                if role=='阳性药物对照':name='RSL3（合成对照）';dose=1
                if key in locations:
                    i,j=locations[key];role='药物库';drug=f'DEMO-{i+1:03}';name=NAMES[i]+'（合成）';dose=5
                layout.append(dict(well=key,role=role,drug_id=drug,name=name,dose=dose,excluded=False,reason='',cell='EV' if cell=='E' else 'H'))
                n0=500*(1+((ord(row)-65)*2+(col==12)-2.5)*.01)
                for t in TIMES:
                    ratio=1
                    if role=='阳性药物对照':ratio=1-.8*t/24
                    if role=='药物库':
                        i,j=locations[key]
                        ratio=CURVES[i][t//2] if cell=='H' else 1-.04*t/24
                        if cell=='E' and i==3:ratio=1-.64*t/24
                        if cell=='H' and i==4:ratio=1-(.85 if j==0 else .15)*t/24
                        ratio*=1+(.015 if j==0 else -.015)*t/24
                    times[str(t)][key]=round(n0*2**(t/24)*ratio,6) if role!='未使用' else 0
        raw=['Vessel Name: SYNTHETIC DEMO 01A'+cell,'Metric: Phase Object Count (Per Image)','Notes: SYNTHETIC DATA ONLY - NOT A REAL EXPERIMENT']
        for t in TIMES:
            raw.extend([f'Time Stamp:\t2026/10/30\tElapsed:\t{t}\thours','\t'+'\t'.join(map(str,range(1,13)))])
            raw.extend(row+'\t'+'\t'.join(str(times[str(t)][f'{row}{c:02}']) for c in range(1,13)) for row in 'ABCDEFGH')
        records.append(dict(plate='01A'+cell,layout=layout,data=dict(header='SYNTHETIC DEMO ONLY',times=times),raw='\n'.join(raw),filename='DEMO_01A'+cell+'.txt',synthetic=True,note='完全合成示例，非真实实验'))
    return records

def write_csv(path,rows):
    cols=list(dict.fromkeys(k for r in rows for k in r))
    with open(path,'w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=cols);writer.writeheader()
        writer.writerows({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in r.items()} for r in rows)

def generate(destination=None):
    import plots,research_plots
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    out=Path(destination or Path(__file__).resolve().parent/'example');out.mkdir(exist_ok=True,parents=True)
    workspace=out/'workspace';workspace.mkdir(exist_ok=True)
    records=synthetic_records()
    store=Store(workspace/'data'/'screen.sqlite')
    existing=store.records()
    expected={r['plate']:(r['raw'],r['layout']) for r in records}
    if existing and {r['plate']:(r['raw'],r['layout']) for r in existing}!=expected:
        store.db.close();raise ValueError('示例工作区已有修改。请使用新的目标目录生成，避免把修改后的数据标记为合成基准。')
    if any(not e.get('synthetic') for e in reviews.history(store)):
        store.db.close();raise ValueError('示例库含额外人工复筛记录。请使用新目标目录生成，不覆盖用户记录。')
    (workspace/'templates.json').write_text(json.dumps({r['plate']:r['layout'] for r in records},ensure_ascii=False),encoding='utf-8')
    if not store.records():
        for r in records:store.save(r)
    records=store.records();profiles=research.temporal_profiles(records);result=research.stability(records,top_k=3)
    result['synthetic']=True
    import qc
    result['source_signature']=json.dumps(dict(records=[(r.get('id'),r['plate'],r.get('hash')) for r in records],baseline=result['baseline'],top_k=3,qc=qc.config()),sort_keys=True)
    old_runs=reviews.runs(store,'stability')
    if not old_runs or old_runs[0]['result'].get('source_signature')!=result['source_signature']:reviews.save_run(store,'stability',result)
    source=[dict(id=r['id'],plate=r['plate'],hash=r['hash']) for r in records]
    if not reviews.latest(store):
        for i,status in [(0,'待复筛'),(1,'待复筛'),(4,'结果不确定')]:
            d=dict(pair='01A',drug_id=f'DEMO-{i+1:03}',name=NAMES[i]+'（合成）',status=status,reviewer='演示记录人',experiment_date='2026-10-30',
                independent_runs=1 if i==4 else 0,experiment_id='SIM-R05' if i==4 else '',evidence='合成记录：模拟复筛仍有分歧，不能确定结果' if i==4 else '',
                notes='仅演示复筛工作流，不代表实际实验',source_records=source,synthetic=True)
            case=reviews.append(store,d)
            if i==0:
                d.update(status='已复现',independent_runs=2,experiment_id='SIM-R01/R02',evidence='合成证据：两次模拟独立实验方向一致；非真实结果',revision_reason='演示从计划到结果的版本追加')
                reviews.append(store,d,case)
    events=reviews.history(store);latest=reviews.latest(store)
    for r in synthetic_records():(out/r['filename']).write_text(r['raw'],encoding='utf-8')
    (out/'synthetic_records.json').write_text(json.dumps(synthetic_records(),ensure_ascii=False,indent=2),encoding='utf-8')
    write_csv(out/'temporal_profiles.csv',profiles);write_csv(out/'stability_summary.csv',result['summary']);write_csv(out/'stability_scenarios.csv',result['details']);write_csv(out/'rescreen_history.csv',events)
    (out/'analysis.json').write_text(json.dumps(dict(synthetic=True,temporal_parameters=research.TEMPORAL_DEFAULTS,profiles=profiles,stability=result,reviews=events),ensure_ascii=False,indent=2),encoding='utf-8')
    fig=plots.Figure(figsize=(13,8.6),dpi=130);FigureCanvasAgg(fig)
    fig.subplots_adjust(left=.065,right=.98,bottom=.09,top=.87,hspace=.68,wspace=.28)
    for i in range(8):
        ax=fig.add_subplot(2,4,i+1);drug=f'DEMO-{i+1:03}'
        for cell,color in [('EV',plots.BLUE),('R273H',plots.PURPLE)]:
            p=next(p for p in profiles if p['drug_id']==drug and p['cell']==cell)
            for rep in p['replicates']:ax.plot(TIMES,rep['curve'],color=color,alpha=.35,ls='--',lw=1)
            ax.plot(TIMES,p['curve'],color=color,label=cell,lw=2)
        hp=next(p for p in profiles if p['drug_id']==drug and p['cell']=='R273H')
        ax.axhline(.8,color='#a77d56',ls=':',lw=1);ax.set_title(f'{drug} {NAMES[i]}\nH：{hp["label"]}',fontsize=10,loc='left')
        ax.set(xlim=(0,24),ylim=(0,1.15),xticks=[0,6,12,18,24],xlabel='Elapsed（h）',ylabel='相对增长 R(t)');ax.grid(alpha=.12)
        if i==0:ax.legend(fontsize=8)
        if i==4:ax.text(.04,.06,'复孔分类不一致',transform=ax.transAxes,color='#b26727',fontsize=9)
    fig.suptitle('合成示例 · 时序反应分型（非真实实验）',fontsize=18,x=.065,ha='left')
    fig.text(.065,.025,'细线为两个技术复孔，粗线为均值；0.8 为示例规则阈值。分类可解释，但不代表机制或临床疗效。',fontsize=10,color='#52666c')
    fig.savefig(out/'temporal_example.png');sf=research_plots.stability_figure(result,8);FigureCanvasAgg(sf);sf.suptitle('合成示例 · 非真实实验',fontsize=13,y=.99);sf.savefig(out/'stability_example.png',dpi=130)
    esc=lambda v:html.escape(str(v))
    hprofiles=[p for p in profiles if p['cell']=='R273H'];stab={(r['pair'],r['drug_id']):r for r in result['summary']}
    table=''.join('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in [p['drug_id'],p['name'],p['label'],p['onset'] if p['onset'] is not None else '未检出',p['agreement'],f"{stab[(p['pair'],p['drug_id'])]['top_hits']}/{len(result['scenarios'])}"] )+'</tr>' for p in sorted(hprofiles,key=lambda p:p['drug_id']))
    rt=''.join('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in [r['case_id'],r['drug_id'],r['status'],r['independent_runs'],r['evidence'] or '尚未复筛'])+'</tr>' for r in latest)
    report=f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>药筛研究模块 · 合成示例</title><style>body{{font:16px/1.65 "Microsoft YaHei",sans-serif;color:#233943;background:#f1f5f6;margin:0}}main{{max-width:1150px;margin:32px auto;background:white;padding:32px;border-radius:16px}}h1,h2{{color:#14566b}}.notice{{padding:16px;background:#fff0d8;border-left:5px solid #d18a2b}}img{{width:100%;height:auto}}table{{width:100%;border-collapse:collapse;font-size:14px}}th,td{{text-align:left;padding:10px;border-bottom:1px solid #dbe3e8}}th{{background:#e9f1f4}}.small{{color:#60727b}}a{{color:#177ca0}}</style><main><h1>三个研究模块的示例结果</h1><p class="notice"><b>全部数据与复筛结论均为合成示例，不是真实药物实验。</b><br>2块96孔板、8种虚构药物、每种细胞两个技术复孔、13个时间点。示例数据库与正式数据库完全分开。</p><h2>1. 时序反应分型</h2><p>先校正各孔0 h和DMSO增长，再按连续低于阈值、出现时间、末端恢复等规则分类。图中保留复孔差异和单点异常。</p><img src="temporal_example.png" alt="合成药物时间响应"><table><tr><th>药物编号</th><th>设计情形</th><th>H均值分类</th><th>规则起点h</th><th>复孔分类</th><th>Top-3入选方案</th></tr>{table}</table><h2>2. 排名稳定性</h2><p>比较{len(result['scenarios'])}个非重复方案，包括指标、综合权重、门槛及单复孔组合。入选频率的分母是全部方案，<b>不是成功概率或统计置信度</b>；灰格保留筛除/不可计算的原因。</p><img src="stability_example.png" alt="各方案排名热图"><h2>3. 复筛记录</h2><p>示例包含计划、结果不确定和追加“已复现”版本。标签均由演示流程人工设定，不是模型预测，也不能用于宣传真实验证成功。</p><table><tr><th>记录</th><th>药物</th><th>当前状态</th><th>模拟独立次数</th><th>模拟证据</th></tr>{rt}</table><h2>检查与复现</h2><p>每条复筛记录保留源数据引用；源板替换后提示复核。稳定性分析另存冻结快照。异常仅提示，不能自动删改结果。</p><p><a href="temporal_profiles.csv">完整时序特征CSV</a> · <a href="stability_summary.csv">稳定性汇总CSV</a> · <a href="stability_scenarios.csv">逐方案CSV</a> · <a href="rescreen_history.csv">复筛版本历史CSV</a> · <a href="analysis.json">完整分析JSON</a></p><p class="small">基础版采用规则和敏感性分析，尚未训练预测模型。后续AI4S需要独立复筛标签、明确基准和不泄漏的验证设计。</p></main></html>'''
    (out/'示例报告.html').write_text(report,encoding='utf-8');store.db.close();return out

if __name__=='__main__':print(generate())
