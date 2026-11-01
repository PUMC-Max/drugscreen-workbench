import plots
import numpy as np
from engine import TIMES

def profile_figure(profile,settings):
    f=plots.Figure(figsize=(10,5.5),dpi=110);ax=f.add_subplot();f.subplots_adjust(left=.1,right=.97,bottom=.17,top=.82)
    for r in profile['replicates']:ax.plot(TIMES,[v if v is not None else np.nan for v in r['curve']],ls='--',lw=1,label=r['well']+' '+r['label'],alpha=.7)
    ax.plot(TIMES,[v if v is not None else np.nan for v in profile['curve']],color=plots.BLUE,lw=2.7,marker='o',ms=3,label='有效复孔均值')
    ax.axhline(settings['inhibition'],color='#b66c30',ls=':',label='抑制阈值')
    ax.axhline(settings['recovery'],color='#6f805a',ls=':',label='恢复阈值')
    if profile['onset'] is not None:ax.axvline(profile['onset'],color='#a0a0a0',lw=1,ls=':')
    ax.set(xlabel='仪器 Elapsed（h）',ylabel='相对增长 R(t)',xticks=TIMES)
    ax.set_title(f"{profile['plate']} · {profile['drug_id']} · {profile['name']}\n{profile['label']} / 复孔：{profile['agreement']}",loc='left',fontsize=12)
    ax.legend(fontsize=8,ncol=2);ax.grid(alpha=.15)
    f.text(.08,.045,'规则分型，不代表机制；'+profile['reason'],fontsize=8,color='#57666e')
    return f

def stability_figure(result,limit=16):
    f=plots.Figure(figsize=(12,6.8),dpi=110);ax=f.add_subplot();f.subplots_adjust(left=.26,right=.94,bottom=.32,top=.87)
    rows=result['summary'][:limit];scenarios=result['scenarios']
    if not rows:
        plots.empty(f,'没有可用于稳定性比较的配对药物');return f
    lookup={(d['pair'],d['drug_id'],d['scenario']):d for d in result['details']}
    matrix=np.array([[lookup[(r['pair'],r['drug_id'],i)]['rank'] or np.nan for i in range(len(scenarios))] for r in rows])
    maxrank=max([d['rank'] for d in result['details'] if d['rank'] is not None]+[2])
    im=ax.imshow(matrix,aspect='auto',cmap=plots.matplotlib.colormaps['viridis_r'].with_extremes(bad='#e7e9eb'),vmin=1,vmax=maxrank)
    ax.set_yticks(range(len(rows)),[f"{r['drug_id']} · {r['name'][:12]}  {r['top_hits']}/{r['total']}" for r in rows],fontsize=8)
    ax.set_xticks(range(len(scenarios)),[s['name'] for s in scenarios],rotation=50,ha='right',fontsize=8)
    for y,r in enumerate(rows):
        for x in range(len(scenarios)):
            d=lookup[(r['pair'],r['drug_id'],x)];value=d['rank']
            rgba=im.cmap(im.norm(value)) if value is not None else (1,1,1,1)
            luminance=.299*rgba[0]+.587*rgba[1]+.114*rgba[2]
            ax.text(x,y,str(value) if value is not None else '筛' if d['state']=='筛除' else '缺',ha='center',va='center',fontsize=7,
                    color='white' if luminance<.48 else '#273a43')
    f.colorbar(im,ax=ax,fraction=.025,pad=.02,label='方案内名次（越小越靠前）')
    ax.set_title(f"Top-{result['top_k']} 排名敏感性 · {len(scenarios)} 个非重复方案",loc='left',fontsize=13)
    f.text(.04,.025,'左侧分数为入选方案数 / 全部方案数，不是成功概率。灰格：筛除或指标不可计算；未修改正式排名。',fontsize=9,color='#57666e')
    return f
