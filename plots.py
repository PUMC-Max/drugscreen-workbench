"""直接使用分析结果绘图；质控仅改变提示样式，不屏蔽数值。"""
import sys, os, math, statistics as st
from pathlib import Path

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'vendor'))
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'data'/'matplotlib'))
import matplotlib
matplotlib.rcParams.update({'font.family':['Microsoft YaHei','SimHei','DejaVu Sans'],
    'axes.unicode_minus':False,'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
    'figure.facecolor':'white','axes.titlepad':14,'svg.fonttype':'path'})
import numpy as np
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle, Patch
from matplotlib.colors import Normalize, TwoSlopeNorm
from engine import analyze, TIMES, ranked_treatment
from plate_model import analysis_records,mixed_cells,cell_name
import qc

SORTS={'综合分':'score','ΔV':'dV','ΔGR':'dGR','H GI':'HGI','H GR':'HGR','H Z':'HZ','H Robust Z':'HRZ','H 时间平均抑制':'HTI'}
BLUE='#237fa0';PURPLE='#8457a0';ORANGE='#d5862e'

def finite(v):return v is not None and math.isfinite(v)

def ordered(rows,metric):
    key=SORTS[metric];lower=key in ('HGR','HZ','HRZ')
    return sorted(rows,key=lambda r:(not finite(r.get(key)),r[key] if lower else -r[key],r['pair'],r['id']) if finite(r.get(key)) else (True,0,r['pair'],r['id']))

def empty(fig,message):
    fig.clear();ax=fig.add_subplot();ax.set_axis_off();ax.text(.5,.5,message,ha='center',va='center',transform=ax.transAxes,color='#60757b',fontsize=13)
    return dict(rows=[],targets={},message=message)

def scatter(fig,rows,metric='V',scope='全部已配对药物',thresholds=None):
    xkey,ykey=('EV','H') if metric=='V' else ('EVGR','HGR')
    valid=[r for r in rows if finite(r.get(xkey)) and finite(r.get(ykey))]
    if not valid:return empty(fig,'暂无可绘制的配对结果\n需要对应 E/H 板以及可计算的 '+metric)
    fig.clear();ax=fig.add_subplot();fig.subplots_adjust(left=.10,right=.95,top=.86,bottom=.15)
    xs=[r[xkey] for r in valid];ys=[r[ykey] for r in valid]
    pts=ax.scatter(xs,ys,s=42,c=BLUE,alpha=.8,edgecolors=[ORANGE if r['flags'] else 'white' for r in valid],linewidths=1.4,picker=6)
    low=min(xs+ys+[0]);high=max(xs+ys+[1]);pad=max((high-low)*.07,.05)
    ax.plot([low-pad,high+pad],[low-pad,high+pad],color='#758891',ls='--',lw=1,label='两组相同')
    ax.set(xlim=(low-pad,high+pad),ylim=(low-pad,high+pad),xlabel='EV · '+metric,ylabel='R273H · '+metric)
    ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.15)
    ax.set_title('EV / R273H 药效对比\n对角线下方：R273H 相对抑制更强',loc='left',fontsize=13)
    if metric=='V' and thresholds:
        if thresholds.get('ev_min') is not None:ax.axvline(thresholds['ev_min'],color='#25805b',ls=':',label='EV 最低 V')
        if thresholds.get('h_max') is not None:ax.axhline(thresholds['h_max'],color='#8d56a3',ls=':',label='H 最高 V')
    ax.legend(loc='upper left',fontsize=8)
    note=f'{scope} · {len(valid)} 点；{len(rows)-len(valid)} 条指标不可计算未绘制。橙框仅表示检查提示。'
    fig.text(.10,.04,note,fontsize=9,color='#52666c')
    data=[dict(pair=r['pair'],drug_id=r['id'],name=r['name'],EV=r[xkey],R273H=r[ykey],metric=metric,flags=r['flags']) for r in valid]
    return dict(rows=data,targets={pts:valid},message=note,ax=ax)

def bars(fig,rows,metric='ΔV',n=15,scope='当前榜单'):
    key=SORTS[metric];valid=[r for r in ordered(rows,metric) if finite(r.get(key))][:n]
    if not valid:return empty(fig,'暂无可绘制的排名结果\n请检查配对数据与当前筛选条件')
    fig.clear();ax=fig.add_subplot();fig.subplots_adjust(left=.37,right=.94,top=.88,bottom=.14)
    y=list(range(len(valid)));values=[r[key] for r in valid]
    rects=ax.barh(y,values,color=[ORANGE if r['flags'] else BLUE for r in valid],height=.65,picker=True)
    names=[f"{i+1:02}  {r['pair']} · {r['name'][:27]}{'…' if len(r['name'])>27 else ''}" for i,r in enumerate(valid)]
    ax.set_yticks(y,names,fontsize=8 if n>20 else 9);ax.invert_yaxis();ax.axvline(0,color='#829297',lw=.7)
    ax.set_xlabel(metric+('（越低排名越靠前）' if key in ('HGR','HZ','HRZ') else '（越高排名越靠前）'))
    ax.set_title(f'{scope} · 前 {len(valid)} 项',loc='left',fontsize=14);ax.grid(axis='x',alpha=.15);ax.set_axisbelow(True)
    ax.margins(x=.18)
    for patch,value in zip(rects,values):
        ax.annotate(f'{value:.3g}',(value,patch.get_y()+patch.get_height()/2),xytext=(4 if value>=0 else -4,0),textcoords='offset points',va='center',ha='left' if value>=0 else 'right',fontsize=8)
    note='点击条目可看时间曲线；橙色为检查提示，数值及排名未被异常检查修改。'
    fig.text(.10,.035,note,fontsize=9,color='#52666c')
    return dict(rows=[dict(rank=i+1,pair=r['pair'],drug_id=r['id'],name=r['name'],metric=metric,value=r[key],flags=r['flags']) for i,r in enumerate(valid)],
                targets={patch:[r] for patch,r in zip(rects,valid)},message=note,ax=ax)

def heat_values(record,metric,time,qc_settings=None):
    report=analyze(record,qc_settings);ts=record['data']['times'];time=str(int(time))
    controls=[w['well'] for w in record['layout'] if w['role']=='DMSO对照' and not w.get('excluded')]
    median=st.median(ts[time][k] for k in controls) if controls else 0
    references={}
    if mixed_cells(record):
        for cell in ('EV','R273H'):
            keys=[w['well'] for w in record['layout'] if cell_name(w)==cell and w['well'] in controls]
            references[cell]=st.median(ts[time][k] for k in keys) if keys else 0
    alerts={k for f in report['qc']['findings'] if f['level']=='需核对' for k in f['wells']}
    rows=[]
    for w in record['layout']:
        k=w['well'];denom=references.get(cell_name(w),0) if references else median
        value=ts[time][k] if metric=='Counts' else ts[time][k]/denom if metric=='V' and denom>0 else None
        if metric not in ('Counts','V'):value=report['metrics'].get(k,{}).get(metric)
        if w['role'] in ('原库空孔','未使用'):value=None
        rows.append(dict(well=k,cell=cell_name(w),role=w['role'],drug_id=w['drug_id'],name=w['name'],metric=metric,time=int(time) if metric in ('Counts','V') else 24,
            value=value,raw_count=ts[time if metric in ('Counts','V') else '24'][k],excluded=w.get('excluded',False),qc_alert=k in alerts,reason=w.get('reason','')))
    return rows

def heatmap(fig,record,metric='V',time=24,qc_settings=None):
    if record is None:return empty(fig,'尚无有效板记录\n保存一块板后即可查看热图，不要求 E/H 配对')
    rows=heat_values(record,metric,time,qc_settings);array=np.full((8,12),np.nan)
    for r in rows:
        if finite(r['value']):array[ord(r['well'][0])-65,int(r['well'][1:])-1]=r['value']
    fig.clear();ax=fig.add_subplot();fig.subplots_adjust(left=.07,right=.89,top=.86,bottom=.16)
    finitevals=array[np.isfinite(array)];low=float(finitevals.min()) if finitevals.size else 0;high=float(finitevals.max()) if finitevals.size else 1
    if low==high:low-=.05;high+=.05
    center=1 if metric=='V' else 0
    norm=TwoSlopeNorm(vmin=low,vcenter=center,vmax=high) if metric!='Counts' and low<center<high else Normalize(low,high)
    cmap=matplotlib.colormaps['viridis' if metric=='Counts' else 'RdBu'].with_extremes(bad='#e7e9eb')
    im=ax.imshow(array,cmap=cmap,norm=norm,aspect='equal');fig.colorbar(im,ax=ax,fraction=.03,pad=.03,label=metric)
    ax.set_xticks(range(12),range(1,13));ax.set_yticks(range(8),list('ABCDEFGH'));ax.xaxis.tick_top()
    ax.set_title(f"{record['plate']} · {metric} · {int(time) if metric in ('Counts','V') else 24} h",loc='left',fontsize=14,pad=28)
    for r in rows:
        x=int(r['well'][1:])-1;y=ord(r['well'][0])-65;v=r['value']
        color='black'
        if finite(v):
            rgba=cmap(norm(v));color='white' if .299*rgba[0]+.587*rgba[1]+.114*rgba[2]<.48 else '#172d36'
        label=f'{v:.3g}' if finite(v) else '空' if r['role'] in ('未使用','原库空孔') else '不可算'
        ax.text(x,y,label,ha='center',va='center',fontsize=8,color=color)
        ax.add_patch(Rectangle((x-.48,y-.48),.96,.96,fill=False,edgecolor='#cf6f1b' if r['qc_alert'] else 'white',lw=2 if r['qc_alert'] else .5))
        if r['excluded']:ax.add_patch(Rectangle((x-.48,y-.48),.96,.96,fill=False,hatch='///',edgecolor='#30383d',lw=.7))
    note='橙框：检查提示（保留数值）；斜线：人工排除（展示原孔值）；灰色：空孔或不可计算。点击孔位看曲线。'
    fig.text(.06,.055,note,fontsize=8,color='#52666c')
    return dict(rows=rows,targets={},message=note,ax=ax,heat={r['well']:r for r in rows},record=record)

def time_series(records,pair,drug_id,mode='V'):
    rows=[]
    for record in analysis_records(records):
        if record['plate'][:-1]!=pair:continue
        group=[w for w in record['layout'] if ranked_treatment(w) and w['drug_id']==drug_id]
        if not group:continue
        ts=record['data']['times'];cell='EV' if record['plate'].endswith('E') else 'R273H'
        controls=[w['well'] for w in record['layout'] if w['role']=='DMSO对照' and not w.get('excluded')]
        for w in group:
            for t in TIMES:
                raw=ts[str(t)][w['well']]
                denom=1 if mode=='Counts' else ts['0'][w['well']] if mode=='相对0h' else st.median(ts[str(t)][k] for k in controls) if controls else 0
                rows.append(dict(plate=record['plate'],cell=cell,well=w['well'],drug_id=drug_id,name=w['name'],time=t,mode=mode,raw_count=raw,value=raw/denom if denom>0 else None,excluded=w.get('excluded',False)))
    return rows

def curves(fig,records,pair,drug_id,mode='V',qc_settings=None):
    records=analysis_records(records)
    rows=time_series(records,pair,drug_id,mode)
    if not rows:return empty(fig,'暂无该药物的有效板记录\n选择一块板中的药物，或点击散点 / 排名条目')
    fig.clear();ax=fig.add_subplot();fig.subplots_adjust(left=.10,right=.96,top=.84,bottom=.17)
    used=[r for r in records if r['plate'][:-1]==pair];points={(r['plate'],k,t) for r in used for f in qc.inspect(r,qc_settings)['findings'] if f['level']=='需核对' for k in f['wells'] for t in f['times']}
    targets={};means=[]
    for cell,color in [('EV',BLUE),('R273H',PURPLE)]:
        wells=sorted({r['well'] for r in rows if r['cell']==cell});group=[r for r in rows if r['cell']==cell]
        for k in wells:
            data=[r for r in group if r['well']==k];ex=data[0]['excluded'];vals=[r['value'] if finite(r['value']) else np.nan for r in data]
            line,=ax.plot(TIMES,vals,color='#969a9c' if ex else color,lw=1,alpha=.7,ls=':' if ex else '--',label=cell+' '+k+('（已排除）' if ex else ''),picker=5)
            targets[line]=data
            flagged=[r for r in data if (r['plate'],k,r['time']) in points and finite(r['value'])]
            ax.scatter([r['time'] for r in flagged],[r['value'] for r in flagged],facecolors='none',edgecolors=ORANGE,s=50,zorder=4)
        vals=[]
        for t in TIMES:
            data=[r for r in group if r['time']==t and not r['excluded']]
            v=st.mean(r['value'] for r in data) if data and all(finite(r['value']) for r in data) else None
            vals.append(v if finite(v) else np.nan);means.append(dict(plate=group[0]['plate'] if group else '',cell=cell,well='有效复孔均值',drug_id=drug_id,name=rows[0]['name'],time=t,mode=mode,value=v,excluded=False))
        if group:ax.plot(TIMES,vals,color=color,lw=2.6,marker='o',markersize=3,label=cell+' 有效复孔均值')
    if mode!='Counts':ax.axhline(1,color='#8b999d',lw=.8,ls=':')
    ax.set(xlabel='仪器 Elapsed（h）',ylabel={'V':'DMSO 归一化 V','相对0h':'N(t) / N(0)','Counts':'每视野平均 counts'}[mode],xticks=TIMES)
    ax.set_title(f'{pair} · {drug_id} · {rows[0]["name"][:55]}\n两种细胞的时间反应',loc='left',fontsize=13)
    ax.grid(alpha=.15);ax.legend(fontsize=8,ncol=2,loc='best')
    cells={r['cell'] for r in rows};note='细线：技术复孔；粗线：有效复孔均值；橙圈：检查提示。没有独立重复，不绘制推断性置信区间。'
    if len(cells)<2:note+=' 另一细胞类型尚缺数据。'
    fig.text(.07,.045,note,fontsize=8,color='#52666c')
    return dict(rows=rows+means,targets=targets,message=note,ax=ax)
