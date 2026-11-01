import tkinter as tk
from tkinter import ttk,filedialog,messagebox
import csv,json
from pathlib import Path
from datetime import datetime
import plots
import overview_plots
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from engine import rankings, VERSION, ranked_treatment
from plate_model import analysis_records,pair_status

class Visualizations:
    def __init__(self,parent,app):
        self.app=app;self.records=[];self.current=[];self.all=[];self.snapshot={};self.drug_map={};self.result=None
        self.kind=tk.StringVar(value='抑制倾向散点');self.scope=tk.StringVar(value='全部已配对药物')
        self.effect=tk.StringVar(value=overview_plots.METRICS[0]);self.page=0
        self.metric=tk.StringVar(value='V');self.sort=tk.StringVar(value='ΔV');self.top=tk.StringVar(value='15')
        self.plate=tk.StringVar();self.heatmetric=tk.StringVar(value='V');self.time=tk.StringVar(value='24')
        self.drug=tk.StringVar();self.curvemode=tk.StringVar(value='V')
        bar=ttk.Frame(parent);bar.pack(fill='x',pady=(0,8))
        self.combo(bar,'图表',self.kind,['抑制倾向散点','抑制差瀑布图','强度与选择性','两组药效散点','选择性排名','96孔热图','单药时间曲线'],18)
        self.scopebox=self.combo(bar,'数据范围',self.scope,['全部已配对药物','当前筛选榜单'],18)
        ttk.Button(bar,text='刷新图表',command=self.render).pack(side='left',padx=8)
        ttk.Button(bar,text='导出图片＋数据',command=self.export).pack(side='right')
        self.controls=ttk.Frame(parent);self.controls.pack(fill='x')
        self.panels={}
        for name in ('抑制倾向散点','抑制差瀑布图','强度与选择性','两组药效散点','选择性排名','96孔热图','单药时间曲线'):self.panels[name]=ttk.Frame(self.controls)
        for name in ('抑制倾向散点','抑制差瀑布图','强度与选择性'):
            self.combo(self.panels[name],'抑制指标',self.effect,overview_plots.METRICS,32)
        for label,delta in [('上一页',-1),('下一页',1)]:ttk.Button(self.panels['抑制差瀑布图'],text=label,command=lambda d=delta:self.turn_page(d)).pack(side='left',padx=4)
        self.pageinfo=ttk.Label(self.panels['抑制差瀑布图'],text='每页30项');self.pageinfo.pack(side='left',padx=8)
        self.combo(self.panels['两组药效散点'],'比较指标',self.metric,['V','GR'])
        self.combo(self.panels['选择性排名'],'排序指标',self.sort,list(plots.SORTS),22)
        self.combo(self.panels['选择性排名'],'显示数量',self.top,['10','15','20','30'],8)
        self.platebox=self.combo(self.panels['96孔热图'],'板号',self.plate,[],30)
        self.combo(self.panels['96孔热图'],'指标',self.heatmetric,['Counts','V','GI','GR','Z','RZ'],12)
        self.timebox=self.combo(self.panels['96孔热图'],'Counts/V 时间 h',self.time,[str(t) for t in plots.TIMES],8)
        self.drugbox=self.combo(self.panels['单药时间曲线'],'药物',self.drug,[],65)
        self.combo(self.panels['单药时间曲线'],'纵轴',self.curvemode,['V','相对0h','Counts'],12)
        self.summary=ttk.Label(parent,text='',wraplength=1160);self.summary.pack(anchor='w',pady=6)
        self.figure=plots.Figure(figsize=(11,6),dpi=100)
        self.canvas=FigureCanvasTkAgg(self.figure,master=parent)
        toolbar_frame=ttk.Frame(parent);toolbar_frame.pack(side='bottom',fill='x')
        self.toolbar=NavigationToolbar2Tk(self.canvas,toolbar_frame,pack_toolbar=False);self.toolbar.pack(side='left')
        self.hover=ttk.Label(parent,text='悬停查看数值；点击散点或条形查看对应药物曲线。',wraplength=1150)
        self.hover.pack(side='bottom',anchor='w',pady=4)
        self.canvas.get_tk_widget().pack(fill='both',expand=True)
        self.canvas.mpl_connect('pick_event',self.pick);self.canvas.mpl_connect('motion_notify_event',self.motion)
        self.canvas.mpl_connect('button_press_event',self.heat_click)
        self.render()
    def combo(self,parent,label,var,values,width=16):
        ttk.Label(parent,text=label).pack(side='left',padx=(2,6))
        box=ttk.Combobox(parent,textvariable=var,values=values,state='readonly',width=width)
        box.pack(side='left',padx=(0,12));box.bind('<<ComboboxSelected>>',lambda e:self.render());return box
    def turn_page(self,delta):
        total=self.result.get('total',0) if self.result else 0
        self.page=max(0,min(self.page+delta,max(0,(total-1)//30)));self.render()
    def update(self,records,current,snapshot):
        self.records=records;self.current=current;self.snapshot=dict(snapshot)
        self.all=rankings(records,weight=float(snapshot.get('weight') or .5),qc_settings=self.app.qc_settings)
        plates=[r['plate'] for r in records];self.platebox.config(values=plates)
        if self.plate.get() not in plates:self.plate.set(plates[0] if plates else '')
        self.drug_map={}
        for r in analysis_records(records):
            for w in r['layout']:
                if ranked_treatment(w):self.drug_map[f"{r['plate'][:-1]} | {w['drug_id']} | {w['name']}"]=(r['plate'][:-1],w['drug_id'])
        names=sorted(self.drug_map);self.drugbox.config(values=names)
        if self.drug.get() not in names:self.drug.set(names[0] if names else '')
        self.render()
    def render(self):
        try:
            kind=self.kind.get()
            for key,panel in self.panels.items():
                panel.pack_forget()
                if key==kind:panel.pack(fill='x',pady=4)
            overview=kind in ('抑制倾向散点','抑制差瀑布图','强度与选择性')
            self.scopebox.config(state='readonly' if overview or kind in ('两组药效散点','选择性排名') else 'disabled')
            self.timebox.config(state='readonly' if self.heatmetric.get() in ('Counts','V') else 'disabled')
            rows=self.current if self.scope.get()=='当前筛选榜单' else self.all
            if kind in ('抑制倾向散点','强度与选择性'):
                self.result=overview_plots.scatter(self.figure,rows,self.effect.get(),self.scope.get(),strength=kind=='强度与选择性')
            elif kind=='抑制差瀑布图':
                self.result=overview_plots.waterfall(self.figure,rows,self.effect.get(),self.scope.get(),offset=self.page*30)
                self.page=self.result.get('offset',0)//30;self.pageinfo.config(text=f"第 {self.page+1} 页 · 共 {self.result.get('total',0)} 项")
            elif kind=='两组药效散点':
                thresholds={k:float(self.snapshot[k]) if self.snapshot.get(k) not in (None,'') else None for k in ('ev_min','h_max')}
                self.result=plots.scatter(self.figure,rows,self.metric.get(),self.scope.get(),thresholds)
            elif kind=='选择性排名':self.result=plots.bars(self.figure,rows,self.sort.get(),int(self.top.get()),self.scope.get())
            elif kind=='96孔热图':
                r=next((r for r in self.records if r['plate']==self.plate.get()),None)
                self.result=plots.heatmap(self.figure,r,self.heatmetric.get(),int(self.time.get()),self.app.qc_settings)
            else:
                pair,drug_id=self.drug_map.get(self.drug.get(),('',''))
                self.result=plots.curves(self.figure,self.records,pair,drug_id,self.curvemode.get(),self.app.qc_settings)
            self.toolbar.update();self.canvas.draw_idle()
            paired,pending=pair_status(self.records)
            self.summary.config(text=f'{len(self.records)} 块有效板 / {paired} 组 EV-R273H 配对 / {len(self.all)} 条可计算配对药效。待配对：'+(', '.join(pending) or '无'))
            self.hover.config(text=self.result['message'])
            self.render_meta=dict(app_version=VERSION,generated_at=datetime.now().astimezone().isoformat(),chart=kind,
                scope=self.scope.get() if overview or kind in ('两组药效散点','选择性排名') else '指定板或药物（不受榜单筛选限制）',
                options=dict(scatter_metric=self.metric.get(),sort=self.sort.get(),top=self.top.get(),plate=self.plate.get(),heat_metric=self.heatmetric.get(),time=self.time.get(),drug=self.drug.get(),curve_mode=self.curvemode.get()),
                effect_metric=self.effect.get(),page=self.page+1,ranking=self.snapshot,qc_parameters=dict(self.app.qc_settings),
                records=[dict(id=r.get('id'),plate=r['plate'],file=r['filename'],hash=r.get('hash'),created=r.get('created')) for r in self.records],
                note=self.result['message'])
        except Exception as e:
            self.result=plots.empty(self.figure,'图表暂未生成：'+str(e));self.canvas.draw_idle();self.render_meta={}
            self.hover.config(text='请检查数据或条件；未修改任何分析结果。')
    def focus(self,pair,drug_id):
        label=next((k for k,v in self.drug_map.items() if v==(pair,drug_id)),None)
        if not label:return
        self.drug.set(label);self.kind.set('单药时间曲线');self.render();self.app.nb.select(self.app.vistab)
    def pick(self,event):
        if self.toolbar.mode or not self.result:return
        rows=self.result.get('targets',{}).get(event.artist)
        if not rows:return
        row=rows[int(event.ind[0]) if hasattr(event,'ind') and len(event.ind) else 0]
        if 'pair' in row:self.focus(row['pair'],row['id'])
    def heat_item(self,event):
        if not self.result or 'heat' not in self.result or event.inaxes!=self.result.get('ax') or event.xdata is None or event.ydata is None:return None
        x=int(math_floor(event.xdata+.5));y=int(math_floor(event.ydata+.5))
        return self.result['heat'].get(f'{chr(65+y)}{x+1:02}') if 0<=x<12 and 0<=y<8 else None
    def heat_click(self,event):
        if self.toolbar.mode or event.button!=1:return
        item=self.heat_item(event)
        if item:self.app.detail(record=self.result['record'],focus=item['well'])
    def motion(self,event):
        item=self.heat_item(event)
        if item:
            v='不可计算 / 空孔' if item['value'] is None else f"{item['value']:.6g}"
            self.hover.config(text=f"{item['well']} | {item['role']} | {item['drug_id']} {item['name']} | {item['metric']}={v} | {'检查提示；' if item['qc_alert'] else ''}{'人工排除：'+item['reason'] if item['excluded'] else ''}");return
        if not self.result or self.toolbar.mode:return
        for artist,rows in self.result.get('targets',{}).items():
            hit,info=artist.contains(event)
            if hit:
                idx=int(info['ind'][0]) if len(info.get('ind',[])) else 0;row=rows[min(idx,len(rows)-1)]
                if 'pair' in row:
                    effects=overview_plots.effects(row,self.effect.get()) if self.kind.get() in ('抑制倾向散点','抑制差瀑布图','强度与选择性') else None
                    detail=f'EV 抑制={effects[0]:.3g}，H 抑制={effects[1]:.3g}，H−EV={effects[1]-effects[0]:+.3g} {overview_plots.unit(self.effect.get())}' if effects else f"EV V={row['EV']:.4g}, H V={row['H']:.4g}, ΔV={row['dV']:.4g}"
                    self.hover.config(text=f"{row['pair']} | {row['id']} | {row['name']} | {detail} | {row['flags'] or '无检查提示'}")
                else:self.hover.config(text=f"{row['cell']} {row['well']} | {row['time']} h | {row['mode']}={row['value']}")
                return
    def export(self):
        if not self.result or not self.result['rows']:return messagebox.showinfo('导出图表','当前没有可导出的图表数据')
        path=filedialog.asksaveasfilename(initialfile='药筛_'+self.kind.get()+'.png',defaultextension='.png',filetypes=[('高清 PNG','*.png'),('矢量 SVG','*.svg')])
        if not path:return
        try:
            target=Path(path);csvpath=target.with_suffix('.data.csv');metapath=target.with_suffix('.meta.json')
            if any(p.exists() for p in (csvpath,metapath)) and not messagebox.askyesno('覆盖导出附表','同名作图数据或参数文件已存在，是否覆盖？'):return
            self.figure.savefig(target,dpi=240,bbox_inches='tight',facecolor='white')
            rows=self.result['rows'];cols=list(dict.fromkeys(k for row in rows for k in row))
            with csvpath.open('w',newline='',encoding='utf-8-sig') as f:
                writer=csv.DictWriter(f,fieldnames=cols);writer.writeheader();writer.writerows(rows)
            metapath.write_text(json.dumps(self.render_meta,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
            self.app.status.set('已导出图片、作图数据 CSV 和参数 JSON')
        except Exception as e:self.app.error(e)

def math_floor(x):
    import math
    return math.floor(x)
