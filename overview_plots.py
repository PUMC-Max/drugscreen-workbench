"""全库抑制倾向与排名；保留正负方向、来源和不可计算计数。"""
import math
import plots
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm

RED='#df5265'
BLUE='#368bbb'
METRICS=['GI（24 h）','计数抑制（24 h）','时间平均抑制（0–24 h）','GR抑制（24 h）']

def effects(row,metric):
    if metric==METRICS[0]:keys=('EVGI','HGI');scale=100;invert=False
    elif metric==METRICS[1]:keys=('EV','H');scale=100;invert=True
    elif metric==METRICS[2]:keys=('EVTI','HTI');scale=100;invert=False
    elif metric==METRICS[3]:keys=('EVGR','HGR');scale=1;invert=True
    else:raise ValueError('未知抑制指标')
    values=[row.get(k) for k in keys]
    if not all(plots.finite(v) for v in values):return None
    return tuple(scale*(1-v if invert else v) for v in values)

def unit(metric):return 'GR单位' if metric==METRICS[3] else '百分点'

def prepare(rows,metric):
    return [dict(row=r,x=v[0],y=v[1],delta=v[1]-v[0]) for r in rows if (v:=effects(r,metric)) is not None]

def frame(fig,title,subtitle):
    fig.clear();fig.set_facecolor('white');ax=fig.add_subplot()
    fig.subplots_adjust(left=.10,right=.88,top=.80,bottom=.17)
    fig.text(.08,.93,title,fontsize=20,weight='bold',color='#172c40')
    fig.text(.08,.87,subtitle,fontsize=10,color='#64758b')
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['bottom','left']].set_color('#d9e2ec');ax.tick_params(colors='#52657a',labelsize=9)
    ax.grid(color='#edf1f5',linewidth=.8);ax.set_axisbelow(True)
    return ax

def footer(fig,text):fig.text(.08,.035,text,fontsize=8,color='#6b7d90')

def scatter(fig,rows,metric=METRICS[0],scope='全部已配对药物',strength=False):
    data=prepare(rows,metric)
    if not data:return plots.empty(fig,'暂无可计算的抑制率配对结果')
    label='抑制强度与选择性' if strength else '药物抑制倾向'
    ax=frame(fig,label,f'{scope}  ·  {len(data)} 项处理  ·  {metric}')
    xs=[d['y'] if strength else d['x'] for d in data];ys=[d['delta'] if strength else d['y'] for d in data]
    bound=max([abs(d['delta']) for d in data]+[.01]);norm=TwoSlopeNorm(vmin=-bound,vcenter=0,vmax=bound)
    cmap=LinearSegmentedColormap.from_list('selectivity',[BLUE,'#f5f4f1',RED])
    pts=ax.scatter(xs,ys,c=[d['delta'] for d in data],cmap=cmap,norm=norm,s=68,edgecolors='#c3cdd8',linewidths=.6,picker=7,zorder=3)
    fig.colorbar(pts,ax=ax,fraction=.035,pad=.025,label='R273H − EV 抑制差（'+unit(metric)+'）')
    ax.axhline(0,color='#a1aebb',lw=.8);ax.axvline(0,color='#a1aebb',lw=.8)
    suffix='（GR单位）' if metric==METRICS[3] else '（%）'
    if strength:
        ax.set_xlabel('R273H 抑制强度'+suffix);ax.set_ylabel('R273H − EV 抑制差（'+unit(metric)+'）')
        ax.margins(.16)
    else:
        low=min(xs+ys+[0]);high=max(xs+ys+[0]);pad=max((high-low)*.16,1 if metric!=METRICS[3] else .05)
        ax.set(xlim=(low-pad,high+pad),ylim=(low-pad,high+pad),xlabel='EV 抑制'+suffix,ylabel='R273H 抑制'+suffix)
        ax.plot([low-pad,high+pad],[low-pad,high+pad],ls=(0,(4,5)),color='#93a3b5',lw=1,zorder=1)
    # 标注优先显示差异最大的处理；拥挤点由悬停和点击定位，不强行堆叠文字。
    occupied=[];shown=0
    for i in sorted(range(len(data)),key=lambda i:-abs(data[i]['delta'])):
        if shown>=12:break
        px,py=ax.transData.transform((xs[i],ys[i]));text=data[i]['row']['name'][:22];width=max(40,len(text)*6)
        for dx,dy in [(9,10),(9,-19),(-width-9,10),(-width-9,-19),(9,27)]:
            box=(px+dx,py+dy,px+dx+width,py+dy+14)
            if any(not(box[2]<b[0] or box[0]>b[2] or box[3]<b[1] or box[1]>b[3]) for b in occupied):continue
            occupied.append(box);ax.annotate(text,(xs[i],ys[i]),xytext=(dx,dy),textcoords='offset pixels',fontsize=8,color='#34495e',
                bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=1));shown+=1;break
    note=f'红：更偏向抑制 R273H；蓝：更偏向抑制 EV。{len(rows)-len(data)} 项指标不可计算；异常提示不屏蔽数值。'
    footer(fig,note+'\n每点为板配对 × 药物；技术复孔取均值，不绘制推断性误差线。悬停查看，点击进入时间曲线。')
    return dict(rows=[dict(pair=d['row']['pair'],drug_id=d['row']['id'],name=d['row']['name'],metric=metric,EV_inhibition=d['x'],H_inhibition=d['y'],difference=d['delta'],flags=d['row']['flags']) for d in data],
                targets={pts:[d['row'] for d in data]},message=note,ax=ax)

def waterfall(fig,rows,metric=METRICS[0],scope='全部已配对药物',offset=0,page_size=30):
    data=sorted(prepare(rows,metric),key=lambda d:(-d['delta'],d['row']['pair'],d['row']['id']))
    if not data:return plots.empty(fig,'暂无可计算的药物排名')
    offset=max(0,min(offset,((len(data)-1)//page_size)*page_size));page=data[offset:offset+page_size]
    ax=frame(fig,'药物抑制瀑布图',f'{scope}  ·  {metric}  ·  排名 {offset+1}–{offset+len(page)} / {len(data)}')
    fig.subplots_adjust(left=.25,right=.93,top=.79,bottom=.17)
    vals=[d['delta'] for d in page];bars=ax.barh(range(len(page)),vals,height=.62,color=[RED if v>0 else BLUE if v<0 else '#adb8c4' for v in vals],picker=True)
    repeated=len({d['row']['name'] for d in data})!=len(data)
    ax.set_yticks(range(len(page)),[f"{offset+i+1:02}   {d['row']['name']}"+(f" · {d['row']['pair']}" if repeated else '') for i,d in enumerate(page)],fontsize=8)
    ax.invert_yaxis();ax.axvline(0,color='#8899aa',lw=.9);ax.grid(axis='y',visible=False)
    low=min([d['delta'] for d in data]+[0]);high=max([d['delta'] for d in data]+[0]);pad=max((high-low)*.14,.02)
    ax.set_xlim(low-pad,high+pad);ax.set_xlabel('R273H 抑制 − EV 抑制（'+unit(metric)+'）')
    for b,v in zip(bars,vals):ax.annotate(f'{v:+.2f}',(v,b.get_y()+b.get_height()/2),xytext=(5 if v>=0 else -5,0),textcoords='offset points',ha='left' if v>=0 else 'right',va='center',fontsize=8,color='#52657a')
    note=f'红：更偏向抑制 R273H；蓝：更偏向抑制 EV。全库排序保留正负结果，{len(rows)-len(data)} 项不可计算。'
    footer(fig,note+'\n全库按板配对 × 药物汇总，不将同名药物跨板强行合并。可翻页查看全部处理。')
    return dict(rows=[dict(rank=offset+i+1,pair=d['row']['pair'],drug_id=d['row']['id'],name=d['row']['name'],metric=metric,value=d['delta'],EV_inhibition=d['x'],H_inhibition=d['y'],flags=d['row']['flags']) for i,d in enumerate(page)],
                targets={b:[d['row']] for b,d in zip(bars,page)},message=note,ax=ax,total=len(data),offset=offset)
