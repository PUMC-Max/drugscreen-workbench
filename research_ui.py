import tkinter as tk
from tkinter import ttk,filedialog,messagebox
import json,csv,copy,threading,queue
from pathlib import Path
import research,reviews,research_plots
from plate_model import analysis_records
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk

class ResearchPanels:
    def __init__(self,app,rootpath):
        self.app=app;self.path=Path(rootpath);self.records=[];self.snapshot={};self.profiles=[];self.result=None;self.busy=False
        from app import table
        self.table=table
        try:self.temporal_settings=research.temporal_config(json.loads((self.path/'temporal_settings.json').read_text(encoding='utf-8')))
        except (OSError,ValueError):self.temporal_settings=research.temporal_config()
        self.temporal_tab=ttk.Frame(app.nb,padding=12);self.stability_tab=ttk.Frame(app.nb,padding=12);self.review_tab=ttk.Frame(app.nb,padding=12)
        for f,name in [(self.temporal_tab,'时序分型'),(self.stability_tab,'排名稳定性'),(self.review_tab,'复筛记录')]:app.nb.add(f,text=name)
        bar=ttk.Frame(self.temporal_tab);bar.pack(fill='x')
        for name,fn in [('设置分型规则',self.temporal_options),('查看分型曲线',self.profile_plot),('导出分型及特征',self.export_profiles),('为所选药物建立复筛记录',self.review_from_profile)]:ttk.Button(bar,text=name,command=fn).pack(side='left',padx=4)
        self.temporal_info=ttk.Label(self.temporal_tab,text='',wraplength=1170);self.temporal_info.pack(anchor='w',pady=9)
        self.temporal_tree=table(self.temporal_tab,['板号','药物编号','名称','类别','起点 h','末3点 R','最低 R','回升幅度','时间抑制','复孔一致性','检查提示'],[80,100,170,120,65,85,85,85,85,110,240])
        self.temporal_tree.bind('<Double-1>',lambda e:self.profile_plot())
        ttk.Label(self.temporal_tab,text='分类是可解释规则，不是训练模型；异常仅提示。灰色或“—”表示未检出起点或特征不可计算，详见判据。').pack(anchor='w',pady=8)
        bar=ttk.Frame(self.stability_tab);bar.pack(fill='x')
        self.metric=tk.StringVar(value='ΔV');self.top=tk.StringVar(value='3' if app.demo else '10')
        for label,var,options in [('基准排序',self.metric,['ΔV','ΔGR','综合分','时间抑制差']),('Top K',self.top,['3','5','10','20','30','50'])]:
            ttk.Label(bar,text=label).pack(side='left',padx=4);box=ttk.Combobox(bar,textvariable=var,values=options,width=14,state='readonly');box.pack(side='left');box.bind('<<ComboboxSelected>>',lambda e:self.mark_stability())
        self.runbutton=ttk.Button(bar,text='计算并保存快照',command=self.start_stability);self.runbutton.pack(side='left',padx=8)
        ttk.Button(bar,text='方案说明',command=lambda:app.text_window('稳定性分析说明',research.STABILITY_HELP)).pack(side='left')
        ttk.Button(bar,text='查看热图',command=self.stability_plot).pack(side='left',padx=5)
        ttk.Button(bar,text='导出本次分析',command=self.export_stability).pack(side='left')
        bar2=ttk.Frame(self.stability_tab);bar2.pack(fill='x',pady=7)
        ttk.Label(bar2,text='历史快照').pack(side='left');self.runvar=tk.StringVar();self.runbox=ttk.Combobox(bar2,textvariable=self.runvar,width=55,state='readonly');self.runbox.pack(side='left',padx=8);self.runbox.bind('<<ComboboxSelected>>',lambda e:self.load_run())
        self.stability_info=ttk.Label(self.stability_tab,text='尚未计算。此分析不会修改正式榜单。',wraplength=1170);self.stability_info.pack(anchor='w',pady=8)
        self.stability_tree=table(self.stability_tab,['配对','药物编号','名称','Top-K 次数','频率','最优名次','中位名次','最差名次','可排序','筛除','不可算'],[70,110,200,100,80,85,85,85,75,75,75])
        self.stability_tree.bind('<Double-1>',lambda e:self.scenario_detail())
        bar=ttk.Frame(self.review_tab);bar.pack(fill='x')
        for name,fn in [('新建复筛记录',self.new_review),('追加记录版本',self.edit_review),('查看版本历史',self.review_history),('导出完整记录',self.export_reviews)]:ttk.Button(bar,text=name,command=fn).pack(side='left',padx=4)
        ttk.Label(self.review_tab,text='结果为人工录入，系统不自动确认实验成功。源数据被替换或移除时提示重新核对；原版本永不覆盖。',wraplength=1170).pack(anchor='w',pady=9)
        self.review_tree=table(self.review_tab,['记录','配对','药物编号','名称','状态','实验编号','日期','独立次数','记录人','源数据状态'],[60,75,110,160,100,120,100,80,100,220])
        self.review_tree.bind('<Double-1>',lambda e:self.review_history())
        self.refresh_runs()
    def update(self,records,snapshot):
        self.records=copy.deepcopy(analysis_records(records));self.snapshot=copy.deepcopy(snapshot)
        self.profiles=research.temporal_profiles(records,self.temporal_settings,self.app.qc_settings)
        self.temporal_tree.delete(*self.temporal_tree.get_children())
        for i,r in enumerate(self.profiles):
            self.temporal_tree.insert('','end',iid=str(i),values=[r['plate'],r['drug_id'],r['name'],r['label']]+[self.f(r[k]) for k in ('onset','end','minimum','rebound','TI')]+[r['agreement'],r['flags']])
        self.temporal_info.config(text=f'{len(self.profiles)} 条细胞类型×药物时序记录。R(t)校正各孔起始值与同板DMSO增长；阈值为探索性规则。')
        self.mark_stability();self.refresh_reviews()
    @staticmethod
    def f(v):return '—' if v is None else f'{v:.3f}' if isinstance(v,float) else v
    def temporal_options(self):
        win=tk.Toplevel(self.app.root);win.title('时序分型规则');win.geometry('650x400');fields={}
        for i,(k,label) in enumerate([('inhibition','持续抑制上限 R'),('recovery','恢复下限 R'),('rebound','最低回升幅度'),('late_hour','迟发起点 h'),('consecutive','连续观测点数')]):
            ttk.Label(win,text=label).grid(row=i,column=0,padx=15,pady=12);v=tk.StringVar(value=str(self.temporal_settings[k]));fields[k]=v;ttk.Entry(win,textvariable=v).grid(row=i,column=1)
        def save():
            try:
                p=research.temporal_config({k:v.get() for k,v in fields.items()});(self.path/'temporal_settings.json').write_text(json.dumps(p),encoding='utf-8')
                self.temporal_settings=p;self.update(self.records,self.snapshot);win.destroy()
            except Exception as e:self.app.error(e)
        ttk.Button(win,text='规则详情',command=lambda:self.app.text_window('时序判据',research.TEMPORAL_HELP)).grid(row=5,column=0,pady=18)
        ttk.Button(win,text='保存并重新分型',command=save).grid(row=5,column=1)
    def selected_profile(self):
        ids=self.temporal_tree.selection()
        if not ids:raise ValueError('请先选择一条时序记录')
        return self.profiles[int(ids[0])]
    def profile_plot(self):
        try:
            p=self.selected_profile();self.show_figure(research_plots.profile_figure(p,self.temporal_settings),'时序分型',dict(parameters=self.temporal_settings,profile=p))
        except Exception as e:self.app.error(e)
    def show_figure(self,figure,title,data):
        win=tk.Toplevel(self.app.root);win.title(title);win.geometry('1200x760')
        def save():
            path=filedialog.asksaveasfilename(initialfile=title+'.png',defaultextension='.png')
            if path:
                try:figure.savefig(path,dpi=220,bbox_inches='tight');Path(path).with_suffix('.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
                except Exception as e:self.app.error(e)
        ttk.Button(win,text='导出图片与完整分析数据',command=save).pack(anchor='e',padx=10,pady=5)
        canvas=FigureCanvasTkAgg(figure,win);toolbar=NavigationToolbar2Tk(canvas,win,pack_toolbar=False);toolbar.pack(side='bottom',fill='x');canvas.get_tk_widget().pack(fill='both',expand=True);canvas.draw()
        win.figure_canvas=canvas
    def export(self,title,rows,payload):
        path=filedialog.asksaveasfilename(initialfile=title+'.csv',defaultextension='.csv')
        if not path:return
        try:
            cols=list(dict.fromkeys(k for r in rows for k in r))
            with open(path,'w',encoding='utf-8-sig',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=cols);writer.writeheader()
                for row in rows:writer.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in row.items()})
            Path(path).with_suffix('.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8');self.app.status.set('已导出CSV与完整JSON')
        except Exception as e:self.app.error(e)
    def export_profiles(self):self.export('时序分型',self.profiles,dict(parameters=self.temporal_settings,profiles=self.profiles,note=research.TEMPORAL_HELP))
    def baseline(self):
        return dict(mode='A' if self.snapshot.get('mode','A').startswith('A') else 'B',metric=self.metric.get(),weight=float(self.snapshot.get('weight') or .5),
                    **{k:float(self.snapshot[k]) if self.snapshot.get(k) not in ('',None) else None for k in ('ev_min','h_max','delta_min')})
    def signature(self):
        return json.dumps(dict(records=[(r.get('id'),r['plate'],r.get('hash')) for r in self.records],baseline=self.baseline(),top_k=int(self.top.get()),qc=self.app.qc_settings),sort_keys=True)
    def mark_stability(self):
        if self.result:
            stale=self.result.get('source_signature')!=self.signature()
            self.stability_info.config(text=('快照已过期：数据、门槛或方案设置已变化，请重算。' if stale else '当前快照与数据、参数一致。')+f" Top-{self.result['top_k']} / {len(self.result['scenarios'])}个方案；频率不是成功概率。双击药物查看逐方案结果。")
    def start_stability(self):
        if self.busy:return
        if not self.records:return self.app.error('请先保存实验记录')
        signature=self.signature();records=copy.deepcopy(self.records);baseline=self.baseline();top=int(self.top.get());settings=dict(self.app.qc_settings)
        self.busy=True;self.runbutton.config(state='disabled');self.stability_info.config(text='正在计算各方案，不修改原始记录或正式排名……');q=queue.Queue()
        def work():
            try:
                result=research.stability(records,baseline,top,settings);result['source_signature']=signature;result['qc_parameters']=settings;q.put((result,None))
            except Exception as e:q.put((None,str(e)))
        threading.Thread(target=work,daemon=True).start()
        def poll():
            try:result,error=q.get_nowait()
            except queue.Empty:self.app.root.after(100,poll);return
            self.busy=False;self.runbutton.config(state='normal')
            if error:return self.app.error(error)
            try:self.accept_result(result)
            except Exception as e:self.app.error(e)
        self.app.root.after(100,poll)
    def accept_result(self,result):
        run_id=reviews.save_run(self.app.store,'stability',result);self.result=result;self.display_stability();self.refresh_runs()
        self.runvar.set(next(k for k,v in self.runmap.items() if v['id']==run_id))
    def display_stability(self):
        self.stability_tree.delete(*self.stability_tree.get_children())
        for i,r in enumerate(self.result['summary']):self.stability_tree.insert('','end',iid=str(i),values=[r['pair'],r['drug_id'],r['name'],f"{r['top_hits']}/{r['total']}",f"{r['frequency']:.0%}"]+[self.f(r[k]) for k in ('rank_min','rank_median','rank_max','eligible','filtered','unavailable')])
        self.mark_stability()
    def refresh_runs(self):
        self.runmap={f"#{r['id']} · {r['created']}":r for r in reviews.runs(self.app.store,'stability')};self.runbox.config(values=list(self.runmap))
        if self.result is None and self.runmap:
            label=next(iter(self.runmap));self.runvar.set(label);self.result=self.runmap[label]['result'];self.display_stability()
    def load_run(self):
        r=self.runmap.get(self.runvar.get())
        if r:self.result=r['result'];self.display_stability()
    def stability_plot(self):
        if self.result:self.show_figure(research_plots.stability_figure(self.result),'排名稳定性',self.result)
    def scenario_detail(self):
        ids=self.stability_tree.selection()
        if not ids:return
        r=self.result['summary'][int(ids[0])];items=[d for d in self.result['details'] if (d['pair'],d['drug_id'])==(r['pair'],r['drug_id'])]
        self.app.text_window(r['drug_id']+' · 逐方案结果','\n'.join(f"{d['scenario_name']}：{d['state']}；名次 {self.f(d['rank'])}；Top-K {'是' if d['top'] else '否'}" for d in items)+'\n\n方案参数：\n'+json.dumps(self.result['scenarios'],ensure_ascii=False,indent=2))
    def export_stability(self):
        if self.result:self.export('排名稳定性',self.result['summary'],self.result)
    def refresh_reviews(self):
        self.review_rows=reviews.latest(self.app.store);self.review_tree.delete(*self.review_tree.get_children())
        for r in self.review_rows:self.review_tree.insert('','end',iid=str(r['case_id']),values=[r['case_id'],r['pair'],r['drug_id'],r['name'],r['status'],r.get('experiment_id',''),r['experiment_date'],r['independent_runs'],r['reviewer'],reviews.source_state(r,self.records)])
    def selected_review(self):
        ids=self.review_tree.selection()
        if not ids:raise ValueError('请先选择复筛记录')
        return next(r for r in self.review_rows if r['case_id']==int(ids[0]))
    def new_review(self):self.review_editor()
    def review_from_profile(self):
        try:p=self.selected_profile();self.review_editor(identity=(p['pair'],p['drug_id']))
        except Exception as e:self.app.error(e)
    def edit_review(self):
        try:self.review_editor(existing=self.selected_review())
        except Exception as e:self.app.error(e)
    def review_history(self):
        try:
            r=self.selected_review();self.app.text_window('复筛记录版本历史',json.dumps(reviews.history(self.app.store,r['case_id']),ensure_ascii=False,indent=2))
        except Exception as e:self.app.error(e)
    def export_reviews(self):self.export('复筛完整记录',reviews.history(self.app.store),dict(events=reviews.history(self.app.store),note='人工复筛标签；无自动验证，包含各版本与源记录引用'))
    def review_editor(self,existing=None,identity=None):
        candidates={f"{p['pair']} | {p['drug_id']} | {p['name']}":(p['pair'],p['drug_id'],p['name']) for p in self.profiles}
        if existing:
            key=f"{existing['pair']} | {existing['drug_id']} | {existing['name']}";candidates[key]=(existing['pair'],existing['drug_id'],existing['name'])
        elif not candidates:return self.app.error('请先保存药物实验数据')
        else:key=next((k for k,v in candidates.items() if v[:2]==identity),next(iter(candidates)))
        win=tk.Toplevel(self.app.root);win.title('追加复筛版本' if existing else '新建复筛记录');win.geometry('850x770')
        selected=tk.StringVar(value=key);ttk.Label(win,text='药物').grid(row=0,column=0,pady=9,padx=12)
        ttk.Combobox(win,textvariable=selected,values=list(candidates),state='disabled' if existing else 'readonly',width=65).grid(row=0,column=1,sticky='w')
        defaults=dict(status='待复筛',experiment_id='',experiment_date=reviews.now()[:10],independent_runs='0',reviewer='',evidence='',evidence_file='',notes='',revision_reason='')
        if existing:defaults.update(existing)
        fields={}
        for i,(field,label) in enumerate([('status','状态'),('experiment_id','复筛实验编号'),('experiment_date','实验日期 YYYY-MM-DD'),('independent_runs','独立实验次数（不含技术复孔）'),('reviewer','记录人'),('evidence','证据 / 结果说明'),('evidence_file','证据文件（可选）'),('notes','观察、处理或复筛计划'),('revision_reason','版本变更原因（更新必填）')],1):
            ttk.Label(win,text=label).grid(row=i,column=0,padx=12,pady=10,sticky='w');var=tk.StringVar(value=str(defaults.get(field,'')));fields[field]=var
            widget=ttk.Combobox(win,textvariable=var,values=reviews.STATUSES,state='readonly',width=56) if field=='status' else ttk.Entry(win,textvariable=var,width=59)
            widget.grid(row=i,column=1,sticky='w')
        def choose():
            p=filedialog.askopenfilename(title='关联本地复筛证据文件')
            if p:fields['evidence_file'].set(p)
        ttk.Button(win,text='选择证据文件',command=choose).grid(row=10,column=0,pady=8)
        rebind=tk.BooleanVar(value=not bool(existing));ttk.Checkbutton(win,text='将本版本绑定到当前有效初筛数据（否则保留原引用）',variable=rebind).grid(row=10,column=1,sticky='w')
        ttk.Label(win,text='状态由人工提交，不会修改初筛排名；已复现/未复现/不确定需有独立实验与证据。\n关联文件会记录路径和SHA256，不上传文件。',wraplength=770).grid(row=11,column=0,columnspan=2,padx=12,pady=15,sticky='w')
        def save():
            try:
                pair,drug,name=candidates[selected.get()];payload={k:v.get().strip() for k,v in fields.items()};payload.update(pair=pair,drug_id=drug,name=name)
                payload['source_records']=[dict(id=r.get('id'),plate=r['plate'],hash=r.get('hash')) for r in self.records if r['plate'][:-1]==pair] if rebind.get() else (existing or {}).get('source_records',[])
                payload['source_features']=[p for p in self.profiles if p['pair']==pair and p['drug_id']==drug] if rebind.get() else (existing or {}).get('source_features',[])
                reviews.append(self.app.store,payload,existing['case_id'] if existing else None);self.refresh_reviews();win.destroy();self.app.nb.select(self.review_tab)
            except Exception as e:self.app.error(e)
        ttk.Button(win,text='保存为新版本',command=save).grid(row=12,column=1,pady=10,sticky='e')
