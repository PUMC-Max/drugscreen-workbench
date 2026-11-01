import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import json, copy, csv, sys
from pathlib import Path
from engine import *
from templates import load_excel

ROOT=Path(__file__).resolve().parent

METHODS='''计算定义与适用范围 · 版本 1.3

固定条件：药物库 5 μM；RSL3 默认 1 μM；每种细胞一次实验、每种药物两个技术复孔。
以仪器 Elapsed=0 为起点。分析的是每视野平均对象计数，不直接代表细胞死亡率。
EV / R273H 分板时按板号前缀和药物编号匹配；同板时按各孔细胞类型分组，不按日期配对。
同板两种细胞各自使用自身DMSO，物理板仍只保存一条记录。日期与原始文件头完整保留。

1  DMSO 归一化
V(t)=本孔计数(t) / 本板同细胞类型有效 DMSO 孔计数(t)的中位数。
默认药效终点为 24 h；V=1 表示与 DMSO 中位数相同。仅细胞孔不混入 DMSO 对照。

2  GI（本系统采用明确的净增长抑制定义）
先计算本孔 f=N24/N0，以及对照 c=各 DMSO 孔 N24/N0 的中位数。
GI=1−(f−1)/(c−1)。GI=0 表示与对照同等增长，GI=1 表示没有净增长。
GI 可以大于 1 或小于 0，不裁剪。它不同于简单的 1−V。
时间辅助指标 TI=1−∫V(t)dt/24，采用 0–24 h 梯形积分；称为时间平均抑制，非浓度 AUC。

3  GR
GR=2^[log2(f)/log2(c)]−1。f=0 时使用极限值 −1。
c≤1 或 N0=0 时不计算 GR 和 GI；0 h 自身的 GR 无定义，本系统只报告 24 h GR。
https://pmc.ncbi.nlm.nih.gov/articles/PMC4887336/

4  Z-score
Z=(本孔 N24−DMSO N24 均值)/DMSO N24 样本标准差。分数越负，计数越低。
参考分布是同板 DMSO 孔，不是该板全部药物。标准差为零时不可计算。

5  Robust Z-score
RZ=(本孔 N24−DMSO N24 中位数)/(1.4826×MAD)。MAD 为中位绝对偏差。
MAD 为零时不可计算，不回退成普通 Z。对照较少时仅作探索性参考。

6  选择性差值
ΔV=EV 的 V−H 的 V；ΔGR=EV 的 GR−H 的 GR。正值表示对 H 相对抑制更强。
先计算每个有效复孔的指标，再取均值；某指标存在无效复孔时该药该指标标记不可计算。
只有对应的 EV/H 细胞分组及相同药物编号都存在时才参与选择性排名。
重复出现在不同模板中的药物按“模板配对＋药物编号”分别报告，不冒充独立重复。

7  配对 t 检验与 BH-FDR：当前设计不启用
两个技术复孔与连续拍摄时间点不是多次独立实验，不输出跨实验推断的 P 值/FDR。
https://www.ncbi.nlm.nih.gov/books/NBK550206/

8  剂量反应与浓度交互模型：不适用
只有 5 μM 一个浓度，不拟合 IC50、GR50 或浓度曲线。显示的是时间曲线。

9  非劣效/等效性检验：当前设计不启用
EV 保留门槛仅作描述性筛选，不声称正式证明等效或非劣。

两类榜单
A：H 抑制且 EV 保留。探索性综合分=100×clip(1−VH)^w×clip(VEV)^(1−w)。
clip 把综合分中的分量限制在 0–1，原始指标保持不变。默认 w=0.5，可调。
A 默认只要求 VH<1、ΔV>0；EV 最低保留门槛由用户另设，未设时不保证 EV 影响小。
B：相对选择性。只要求 ΔV>0，默认按 ΔV 降序，允许 EV 同时受影响。
两榜可以重叠；可切换按 ΔGR、H 的 GI/GR/Z/RZ/TI 等单项指标排序。
筛选门槛留空表示不限制；所有门槛、综合分权重都是排序规则，不是已验证的生物学标准。
结果含质控提示，不自动丢弃离群孔。手动排除必须填写原因并保留原始数据。

记录与恢复
全部数据存放于 data/screen.sqlite；每次新增、替换、移除前自动生成 data/backups 备份。
移除操作只撤出当前分析，历史记录仍可查看。上传新文件替换会保留旧版本。
备份数据库包含原始 TXT、板位快照和参数之外的原始记录，可复制到其他位置保管。
当前排名参数保存在 settings.json；导出的 CSV 会记录当次筛选参数。
'''

def table(parent,columns,widths=None):
    frame=ttk.Frame(parent); frame.pack(fill='both',expand=True)
    tree=ttk.Treeview(frame,columns=columns,show='headings',selectmode='browse')
    for i,c in enumerate(columns):
        tree.heading(c,text=c); tree.column(c,width=(widths[i] if widths else 100),minwidth=55,anchor='w')
    vs=ttk.Scrollbar(frame,orient='vertical',command=tree.yview)
    hs=ttk.Scrollbar(frame,orient='horizontal',command=tree.xview)
    tree.configure(yscrollcommand=vs.set,xscrollcommand=hs.set)
    tree.grid(row=0,column=0,sticky='nsew');vs.grid(row=0,column=1,sticky='ns');hs.grid(row=1,column=0,sticky='ew')
    frame.rowconfigure(0,weight=1);frame.columnconfigure(0,weight=1)
    return tree

def fmt(v):return '不可计算' if v is None else f'{v:.4f}' if isinstance(v,float) else str(v)

class App:
    def __init__(self,root,demo=False,workspace=None,example_label=None):
        self.demo=demo
        self.workspace=Path(workspace or ROOT)
        self.example_label=example_label or ('合成示例，非真实实验' if demo else '')
        self.root=root;root.title((self.example_label+' · ' if self.example_label else '')+'药筛工作台 · EV / R273H');root.geometry('1320x880');root.minsize(1050,720)
        import theme
        theme.apply(root)
        root.geometry('1500x950');root.minsize(1250,780)
        self.templates=json.loads((self.workspace/'templates.json').read_text(encoding='utf-8'))
        self.store=Store(self.workspace/'data'/'screen.sqlite');self.loaded=None;self.layout=[]
        try:self.qc_settings=qc.config(json.loads((self.workspace/'qc_settings.json').read_text(encoding='utf-8')))
        except (FileNotFoundError,ValueError):self.qc_settings=qc.config()
        ttk.Label(root,text='药筛工作台'+(' · '+self.example_label if self.example_label else ''),style='Title.TLabel').pack(anchor='w',padx=22,pady=(18,3))
        ttk.Label(root,text='5 μM 单浓度筛选  /  0–24 h 连续计数  /  本地保存',foreground='#58717b').pack(anchor='w',padx=24,pady=(0,15))
        body=ttk.Frame(root);body.pack(fill='both',expand=True,padx=18,pady=8)
        sidebar=tk.Frame(body,bg='#172c40',width=195);sidebar.pack(side='left',fill='y',padx=(0,18));sidebar.pack_propagate(False)
        self.nb=ttk.Notebook(body,style='Workspace.TNotebook');self.nb.pack(side='left',fill='both',expand=True)
        self.newtab=ttk.Frame(self.nb,padding=12);self.histtab=ttk.Frame(self.nb,padding=12);self.ranktab=ttk.Frame(self.nb,padding=12);self.methtab=ttk.Frame(self.nb,padding=12)
        for f,t in [(self.newtab,'新增 / 替换读板'),(self.histtab,'实验记录'),(self.ranktab,'两类药物榜单'),(self.methtab,'方法与适用范围')]:self.nb.add(f,text=t)
        self.qctab=ttk.Frame(self.nb,padding=12);self.nb.insert(2,self.qctab,text='异常检查')
        self.new_ui();self.history_ui();self.rank_ui()
        self.qc_ui()
        self.vistab=ttk.Frame(self.nb,padding=12);self.nb.insert(4,self.vistab,text='结果可视化')
        from visualization import Visualizations
        self.visual=Visualizations(self.vistab,self)
        from research_ui import ResearchPanels
        self.research=ResearchPanels(self,self.workspace)
        theme.navigation(sidebar,self.nb)
        tx=tk.Text(self.methtab,wrap='word',font=('Microsoft YaHei UI',11),padx=18,pady=12);tx.pack(fill='both',expand=True);tx.insert('1.0',METHODS);tx.config(state='disabled')
        self.status=tk.StringVar(value=self.example_label+'工作区：与正式数据分开保存，使用完整工作台。' if self.example_label else '本地工作区：数据保存在本机，异常仅提示。')
        ttk.Label(root,textvariable=self.status).pack(anchor='w',padx=22,pady=8)
        self.refresh()
        self.nb.select(self.vistab)
        root.protocol('WM_DELETE_WINDOW',self.close)
    def close(self):
        self.store.db.close();self.root.destroy()
    def error(self,e):messagebox.showerror('未完成操作',str(e),parent=self.root)
    def new_ui(self):
        bar=ttk.Frame(self.newtab);bar.pack(fill='x')
        ttk.Label(bar,text='① 选择板号').pack(side='left')
        self.plate=tk.StringVar(value=next(iter(self.templates)))
        cb=ttk.Combobox(bar,textvariable=self.plate,values=list(self.templates),state='readonly',width=29);cb.pack(side='left',padx=9);cb.bind('<<ComboboxSelected>>',lambda e:self.load_template())
        ttk.Button(bar,text='② 选择 TXT 读板文件',command=self.import_file).pack(side='left',padx=10)
        ttk.Button(bar,text='查看文件头 / 时间点',command=self.preview).pack(side='left')
        ttk.Button(bar,text='更新默认板位 Excel',command=self.update_templates).pack(side='right')
        self.platebox=cb
        self.filelabel=ttk.Label(self.newtab,text='未选择读板文件');self.filelabel.pack(anchor='w',pady=8)
        qb=ttk.Frame(self.newtab);qb.pack(fill='x')
        self.draftqc=ttk.Label(qb,text='导入后自动检查；红色边框表示有待核对提示。');self.draftqc.pack(side='left')
        ttk.Button(qb,text='查看本次异常清单',command=self.draft_report).pack(side='right')
        ttk.Label(self.newtab,text='③ 核对实际板位：点击孔位或双击下表修改细胞、药物和对照信息。修改仅影响本次记录。').pack(anchor='w')
        self.gridframe=ttk.Frame(self.newtab);self.gridframe.pack(fill='x',pady=9)
        self.wellbuttons={}
        for i,r in enumerate('ABCDEFGH'):
            for c in range(1,13):
                key=f'{r}{c:02}'
                b=tk.Button(self.gridframe,text=key,command=lambda k=key:self.edit_well(k),font=('Microsoft YaHei UI',9),relief='flat')
                b.grid(row=i,column=c-1,padx=2,pady=2,sticky='ew');self.wellbuttons[key]=b
        for c in range(12):self.gridframe.columnconfigure(c,weight=1)
        self.welltree=table(self.newtab,['孔位','细胞','类型','药物编号','药物名称','浓度 μM','排除','0 h','24 h'],[65,75,110,100,240,85,65,90,90])
        self.welltree.bind('<Double-1>',lambda e:self.edit_well(self.welltree.selection()[0]) if self.welltree.selection() else None)
        bottom=ttk.Frame(self.newtab);bottom.pack(fill='x',pady=10)
        ttk.Label(bottom,text='实验备注').pack(side='left');self.note=tk.StringVar();ttk.Entry(bottom,textvariable=self.note,width=45).pack(side='left',padx=8)
        self.confirm=tk.BooleanVar();ttk.Checkbutton(bottom,text='已核对板号、文件和实际加药',variable=self.confirm).pack(side='left',padx=8)
        ttk.Button(bottom,text='④ 保存并重新计算',command=self.save).pack(side='right')
        self.load_template()
    def update_templates(self):
        path=filedialog.askopenfilename(title='选择新版默认板位表',filetypes=[('Excel','*.xlsx')])
        if not path:return
        try:
            templates=load_excel(path,default_plate=self.plate.get())
            if not messagebox.askyesno('更新默认模板',f'读取到 {len(templates)} 块板。更新默认布局并清空当前未保存的导入？历史记录保持各自已保存的布局。'):return
            self.templates=templates
            (self.workspace/'templates.json').write_text(json.dumps(templates,ensure_ascii=False),encoding='utf-8')
            self.platebox.config(values=list(templates))
            if self.plate.get() not in templates:self.plate.set(next(iter(templates)))
            self.load_template();self.status.set('默认板位模板已更新，历史记录未改变')
        except Exception as e:self.error(e)
    def load_template(self):
        self.layout=copy.deepcopy(self.templates[self.plate.get()]);self.confirm.set(False);self.loaded=None
        self.filelabel.config(text='未选择读板文件');self.render_wells()
    def render_wells(self):
        self.welltree.delete(*self.welltree.get_children())
        colors={'药物库':'#e7f1f5','DMSO对照':'#cee9df','阳性药物对照':'#ffe0cb','仅细胞':'#e8e2fa','参考药物':'#e1e8fb','原库空孔':'#eeeeee','未使用':'#eeeeee'}
        flagged=set()
        if self.loaded:
            report=qc.inspect(self.draft_record(),self.qc_settings)
            flagged={k for f in report['findings'] if f['level']=='需核对' for k in f['wells']}
            self.draftqc.config(text=f"本次检查：{report['alerts']} 条需核对，{report['unavailable']} 项无法评估；不会自动排除数据。")
        else:self.draftqc.config(text='导入后自动检查；红色边框表示有待核对提示。')
        for w in self.layout:
            k=w['well'];a=b='—'
            if self.loaded:a=self.loaded['data']['times']['0'][k];b=self.loaded['data']['times']['24'][k]
            self.welltree.insert('', 'end',iid=k,values=(k,cell_name(w) or '—',w['role'],w['drug_id'],w['name'],w['dose'] if w['dose'] is not None else '未提供','是' if w.get('excluded') else '',fmt(a),fmt(b)))
            self.wellbuttons[k].config(bg='#d0d0d0' if w.get('excluded') else colors[w['role']],text=k+' '+({'DMSO对照':'DMSO','阳性药物对照':'阳性','仅细胞':'细胞','药物库':'药物','参考药物':'参考'}.get(w['role'],'空')))
            self.wellbuttons[k].config(highlightthickness=2,highlightbackground='#b94336' if k in flagged else '#eeeeee')
    def edit_well(self,key):
        w=next(w for w in self.layout if w['well']==key)
        win=tk.Toplevel(self.root);win.title(f'实际板位 · {key}');win.geometry('550x510');win.transient(self.root)
        fields={}
        for row,(label,k) in enumerate([('细胞类型','cell'),('孔类型','role'),('药物编号','drug_id'),('药物名称','name'),('浓度 μM（可留空）','dose'),('排除原因 / 备注','reason')]):
            ttk.Label(win,text=label).grid(row=row,column=0,padx=15,pady=12,sticky='w')
            value=(cell_name(w) or ('EV' if self.plate.get().endswith('E') else 'R273H')) if k=='cell' else w.get(k,'')
            v=tk.StringVar(value='' if value is None else str(value));fields[k]=v
            ent=ttk.Combobox(win,textvariable=v,values=ROLES if k=='role' else ['EV','R273H'],state='readonly') if k in ('role','cell') else ttk.Entry(win,textvariable=v,width=35)
            ent.grid(row=row,column=1,padx=8)
        ranked=tk.BooleanVar(value=bool(ranked_treatment(w)));ttk.Checkbutton(win,text='纳入药物榜单（阳性/参考药物可选）',variable=ranked).grid(row=6,column=1,pady=8)
        ex=tk.BooleanVar(value=w.get('excluded',False));ttk.Checkbutton(win,text='排除此孔（保留原始数据）',variable=ex).grid(row=7,column=1,pady=8)
        def commit():
            try:
                updated={k:v.get().strip() for k,v in fields.items()};updated['dose']=float(updated['dose']) if updated['dose'] else None
                updated['include_in_ranking']=ranked.get()
                if updated['dose'] is not None and (not math.isfinite(updated['dose']) or updated['dose']<0):raise ValueError('浓度必须为有限非负数')
                if ex.get() and not updated['reason']:raise ValueError('排除孔必须说明原因')
                candidate=copy.deepcopy(self.layout);idx=next(i for i,v in enumerate(candidate) if v['well']==key)
                candidate[idx].update(updated,excluded=ex.get());validate_layout(candidate)
                self.layout=candidate;self.confirm.set(False);self.render_wells();win.destroy()
            except Exception as e:self.error(e)
        ttk.Button(win,text='保存本孔修改',command=commit).grid(row=8,column=1,pady=15)
    def import_file(self):
        path=filedialog.askopenfilename(title='选择整板 0–24 h counts 文件',filetypes=[('计数文件','*.txt')])
        if not path:return
        try:
            raw,data=read_txt(path);self.loaded=dict(raw=raw,data=data,filename=Path(path).name)
            self.filelabel.config(text=f'{Path(path).name}  ·  13 时间点 / 96 孔 / 1,248 数值  ·  请核对文件头与板号')
            self.confirm.set(False);self.render_wells();self.preview()
        except Exception as e:self.error(e)
    def preview(self):
        if not self.loaded:return messagebox.showinfo('文件预览','请先选择 TXT 文件')
        self.text_window('文件信息',self.loaded['data']['header']+'\n\n时间点：'+', '.join(self.loaded['data']['times'])+' h\n\n文件头日期不用于 E/H 配对。请人工确认该文件属于所选板号。')
    def text_window(self,title,content):
        win=tk.Toplevel(self.root);win.title(title);win.geometry('960x650');tx=tk.Text(win,wrap='word',font=('Microsoft YaHei UI',11),padx=16,pady=12);tx.pack(fill='both',expand=True);tx.insert('1.0',content);tx.config(state='disabled')
    def save(self):
        try:
            if not self.loaded:raise ValueError('请导入 TXT 文件')
            if not self.confirm.get():raise ValueError('请先核对并勾选确认')
            record=self.draft_record()
            validate_layout(self.layout);report=analyze(record,self.qc_settings)
            warnings=list(dict.fromkeys(report['warnings']+qc.summaries(report['qc'])))
            if warnings and not messagebox.askyesno('质控提示','\n'.join(warnings[:12])+f"\n共 {report['qc']['alerts']} 条需核对。完整定位请查看本次异常清单。\n\n仍然保存并保留提示？"):return
            existing=any(r['plate']==record['plate'] for r in self.store.records())
            if existing and not messagebox.askyesno('替换读板结果','该板已有有效数据。替换后旧版本保留在历史中，仅新版本参加计算。继续？'):return
            self.store.save(record,replace=existing);self.refresh();self.confirm.set(False)
            self.status.set(f"已保存 {record['plate']}；历史与榜单已更新。");self.nb.select(self.histtab)
        except Exception as e:self.error(e)
    def history_ui(self):
        bar=ttk.Frame(self.histtab);bar.pack(fill='x',pady=(0,10))
        for label,cmd in [('查看计数矩阵 / 单孔曲线',self.detail),('查看单板药效',self.plate_results),('移除所选读板结果',self.remove),('导出所选原始 TXT',self.export_raw),('备份全部记录',self.export_backup)]:ttk.Button(bar,text=label,command=cmd).pack(side='left',padx=3)
        self.history=table(self.histtab,['记录','板号','状态','保存时间','原文件','备注'],[65,90,110,170,300,300])
        self.history.bind('<Double-1>',lambda e:self.detail())
        extra=ttk.Frame(self.histtab);extra.pack(fill='x',pady=5)
        ttk.Button(extra,text='查看所选记录异常',command=self.history_report).pack(side='left')
        ttk.Button(extra,text='载入所选记录核对 / 排除孔',command=self.load_history).pack(side='left',padx=8)
        ttk.Label(self.histtab,text='移除会撤出当前分析并自动更新排名，旧记录仍保留。可重新选板上传新 TXT。').pack(anchor='w',pady=10)
    def selected(self):
        ids=self.history.selection()
        if not ids:raise ValueError('请先选择一条实验记录')
        return next(r for r in self.store.records(False) if str(r['id'])==ids[0])
    def remove(self):
        try:
            r=self.selected()
            if not r['active']:raise ValueError('该记录已不参与分析')
            if messagebox.askyesno('移除读板结果',f"将 {r['plate']} 的记录 {r['id']} 移出分析？对应选择性结果也会更新。"):
                self.store.remove(r['id']);self.refresh();self.status.set('已移除读板结果并重新计算')
        except Exception as e:self.error(e)
    def export_raw(self):
        try:
            r=self.selected();p=filedialog.asksaveasfilename(initialfile=r['filename'],defaultextension='.txt')
            if p:Path(p).write_text(r['raw'],encoding='utf-8')
        except Exception as e:self.error(e)
    def export_backup(self):
        p=filedialog.asksaveasfilename(initialfile='药筛全部记录.sqlite',defaultextension='.sqlite')
        if p:
            if Path(p).resolve()==self.store.path.resolve():return self.error('请选择当前数据库以外的备份位置')
            dest=sqlite3.connect(p);self.store.db.backup(dest);dest.close();self.status.set('完整数据库备份已导出')
    def detail(self,record=None,focus=None,highlight_times=None):
        try:r=record if record is not None else self.selected()
        except Exception as e:return self.error(e)
        win=tk.Toplevel(self.root);win.title(r['plate']+' · 历史计数与曲线');win.geometry('1160x800')
        ttk.Label(win,text=f"记录 {r.get('id','未保存')} / {r['filename']} / {'有效' if r.get('active') else '草稿或历史版本'}；红点为选定提示的时间点").pack(anchor='w',padx=12,pady=8)
        tr=table(win,['孔位','类型','药物','排除']+[f'{t} h' for t in TIMES],[60,100,180,50]+[85]*13)
        canvas=tk.Canvas(win,height=250,bg='white');canvas.pack(fill='x',padx=12,pady=12)
        layout={w['well']:w for w in r['layout']}
        for k,w in layout.items():tr.insert('','end',iid=k,values=[k,w['role'],w['name'],w.get('reason','')]+[fmt(r['data']['times'][str(t)][k]) for t in TIMES])
        def draw(e=None):
            if not tr.selection():return
            key=tr.selection()[0];values=[r['data']['times'][str(t)][key] for t in TIMES]
            c=[w['well'] for w in r['layout'] if w['role']=='DMSO对照' and not w.get('excluded') and (not mixed_cells(r) or cell_name(w)==cell_name(layout[key]))]
            controls=[st.median(r['data']['times'][str(t)][w] for w in c) for t in TIMES] if c else []
            peers=[w['well'] for w in r['layout'] if ranked_treatment(w) and ranked_treatment(layout[key]) and w['drug_id']==layout[key]['drug_id'] and w['well']!=key and not w.get('excluded') and (not mixed_cells(r) or cell_name(w)==cell_name(layout[key]))]
            paired=[r['data']['times'][str(t)][peers[0]] for t in TIMES] if peers else []
            canvas.delete('all');width=max(canvas.winfo_width(),800);maxv=max(values+controls+paired+[1])*1.12
            canvas.create_text(65,15,anchor='w',text=f'{key} · {layout[key]["name"]}   蓝：本孔   绿：DMSO 中位数   紫：另一复孔 {peers[0] if peers else "无"}（原始 counts）')
            canvas.create_line(65,40,65,205,width-25,205,fill='#87979e')
            for n in range(5):
                y=205-n*40;canvas.create_text(57,y,anchor='e',text=f'{maxv*n/4:.0f}')
            for series,color in [(values,'#237aab'),(controls,'#339777'),(paired,'#9461aa')]:
                if not series:continue
                pts=[]
                for t,v in zip(TIMES,series):pts.extend([65+t/24*(width-100),205-v/maxv*160])
                canvas.create_line(*pts,fill=color,width=2)
            for t in (highlight_times or []):
                if t in TIMES:
                    x=65+t/24*(width-100);y=205-values[TIMES.index(t)]/maxv*160
                    canvas.create_oval(x-4,y-4,x+4,y+4,fill='#c13530',outline='')
            for t in TIMES:canvas.create_text(65+t/24*(width-100),225,text=str(t))
            canvas.create_text(width-40,245,text='小时')
        tr.bind('<<TreeviewSelect>>',draw);canvas.bind('<Configure>',draw);tr.selection_set(focus or 'A01');tr.see(focus or 'A01')
    def plate_results(self):
        try:r=self.selected();a=analyze(r,self.qc_settings)
        except Exception as e:return self.error(e)
        win=tk.Toplevel(self.root);win.title(r['plate']+' · 单板描述性结果');win.geometry('1200x650')
        ttk.Label(win,text='；'.join(a['warnings']) or '已计算。该结果不含跨实验显著性推断。',wraplength=1100).pack(anchor='w',padx=12,pady=12)
        cols=['细胞','药物编号','名称','孔位','V','GI','GR','Z','Robust Z','时间平均抑制','提示'];tr=table(win,cols,[85,100,230,100]+[85]*6+[250])
        for d in a['drugs']:tr.insert('','end',values=[d.get('cell',cell_name(r['layout'][0]) or '—')]+[fmt(d[k]) for k in ['id','name','wells','V','GI','GR','Z','RZ','TI','flags']])
    def rank_ui(self):
        bar=ttk.Frame(self.ranktab);bar.pack(fill='x')
        self.mode=tk.StringVar(value='A · H 抑制与 EV 保留');self.sort=tk.StringVar(value='综合分')
        ttk.Combobox(bar,textvariable=self.mode,values=['A · H 抑制与 EV 保留','B · 相对选择性'],state='readonly',width=24).pack(side='left',padx=4)
        ttk.Combobox(bar,textvariable=self.sort,values=['综合分','ΔV','ΔGR','H GI','H GR','H Z','H Robust Z','H 时间平均抑制'],state='readonly',width=20).pack(side='left',padx=4)
        ttk.Button(bar,text='应用条件 / 更新排名',command=self.refresh_ranks).pack(side='left',padx=5)
        ttk.Button(bar,text='导出当前榜单 CSV',command=self.export_ranks).pack(side='right')
        filters=ttk.Frame(self.ranktab);filters.pack(fill='x',pady=12)
        self.params={}
        saved={}
        try:saved=json.loads((self.workspace/'settings.json').read_text(encoding='utf-8'))
        except (FileNotFoundError,json.JSONDecodeError):pass
        for label,key,default in [('EV 最低 V','ev_min',''),('H 最高 V','h_max',''),('最低 ΔV','delta_min',''),('综合分 H 权重 0–1','weight','0.5')]:
            ttk.Label(filters,text=label).pack(side='left',padx=5);v=tk.StringVar(value=saved.get(key,default));self.params[key]=v;ttk.Entry(filters,textvariable=v,width=8).pack(side='left')
        self.rankinfo=ttk.Label(self.ranktab,text='',wraplength=1180);self.rankinfo.pack(anchor='w',pady=(0,10))
        self.rankcols=['排名','配对','药物编号','名称','综合分','EV V','H V','ΔV','ΔGR','H GI','H GR','H Z','H Robust Z','H 时间平均抑制','提示']
        self.ranktree=table(self.ranktab,self.rankcols,[55,65,100,230]+[95]*10+[260])
        self.ranktree.bind('<Double-1>',lambda e:self.chart_selected())
        self.mode.trace_add('write',lambda *a:self.change_mode());self.sort.trace_add('write',lambda *a:self.refresh_ranks())
    def change_mode(self):
        self.sort.set('综合分' if self.mode.get().startswith('A') else 'ΔV')
    def chart_selected(self):
        ids=self.ranktree.selection()
        if ids and hasattr(self,'visual'):
            r=self.rank_data[self.ranktree.index(ids[0])];self.visual.focus(r['pair'],r['id'])
    def draft_record(self):
        if not self.loaded:raise ValueError('请先导入 TXT 文件')
        return dict(self.loaded,plate=self.plate.get(),layout=self.layout,note=self.note.get(),qc_settings=dict(self.qc_settings))
    def draft_report(self):
        try:self.show_qc(self.draft_record(),draft=True)
        except Exception as e:self.error(e)
    def history_report(self):
        try:self.show_qc(self.selected())
        except Exception as e:self.error(e)
    def load_history(self):
        try:
            r=self.selected()
            if self.loaded and not messagebox.askyesno('载入历史记录','将替换当前未保存的导入和板位修改，继续？'):return
            self.plate.set(r['plate']);self.layout=copy.deepcopy(r['layout'])
            self.loaded={k:r[k] for k in ('raw','data','filename','provenance','example','synthetic') if k in r};self.note.set(r.get('note',''))
            self.confirm.set(False);self.filelabel.config(text=r['filename']+' · 从历史记录载入，请核对后保存为新版本')
            self.render_wells();self.nb.select(self.newtab)
        except Exception as e:self.error(e)
    def qc_ui(self):
        top=ttk.Frame(self.qctab);top.pack(fill='x',pady=(0,8))
        ttk.Button(top,text='检查阈值 / 规则说明',command=self.qc_options).pack(side='left')
        ttk.Button(top,text='按当前阈值重查',command=self.refresh_qc).pack(side='left',padx=8)
        ttk.Button(top,text='查看所选异常与曲线',command=self.qc_open_selected).pack(side='left')
        ttk.Button(top,text='导出当前异常清单',command=self.qc_export_all).pack(side='right')
        self.qcfilter=tk.StringVar(value='全部');ttk.Combobox(top,textvariable=self.qcfilter,values=['全部','需核对','无法评估','信息','已排除'],state='readonly',width=12).pack(side='right',padx=8)
        self.qcfilter.trace_add('write',lambda *a:self.refresh_qc())
        self.qcsummary=ttk.Label(self.qctab,text='',wraplength=1170);self.qcsummary.pack(anchor='w',pady=8)
        self.qccols=['记录','板号','级别','检查项','孔位','时间 h','原因']
        self.qctree=table(self.qctab,self.qccols,[60,80,85,170,160,100,600])
        self.qctree.tag_configure('需核对',foreground='#a5372c');self.qctree.tag_configure('无法评估',foreground='#8b650f')
        self.qctree.bind('<Double-1>',lambda e:self.qc_open_selected())
        ttk.Label(self.qctab,text='仅提示：触发阈值不会把结果变为 NA，不改值、不排除、不改变排名。旧版本在“实验记录”查看。').pack(anchor='w',pady=10)
    def refresh_qc(self):
        if not hasattr(self,'qctree'):return
        self.qctree.delete(*self.qctree.get_children());self.qcitems={};self.qcrows=[]
        records=self.store.records();alerts=unknown=0
        for r in records:
            report=qc.inspect(r,self.qc_settings);alerts+=report['alerts'];unknown+=report['unavailable']
            for f in report['findings']:
                if self.qcfilter.get()!='全部' and f['level']!=self.qcfilter.get():continue
                row=[r['id'],r['plate'],f['level'],f['title'],','.join(f['wells']),','.join(map(str,f['times'])),f['reason']]
                iid=self.qctree.insert('','end',values=row,tags=(f['level'],));self.qcitems[iid]=(r,f);self.qcrows.append(row)
        self.qcsummary.config(text=f'{len(records)} 块有效板：{alerts} 条需核对，{unknown} 项无法评估。阈值为探索性设置；没有提示也不等于实验合格。')
    def qc_open_selected(self):
        ids=self.qctree.selection()
        if not ids:return self.error('请先选择一条异常或检查记录')
        r,f=self.qcitems[ids[0]];self.show_qc(r)
        if f['wells']:self.detail(record=r,focus=f['wells'][0],highlight_times=f['times'])
    def qc_export_all(self):
        self.export_qc(self.qccols,self.qcrows,{'parameters':self.qc_settings,'scope':'当前有效记录','filter':self.qcfilter.get(),'qc_version':qc.QC_VERSION})
    def export_qc(self,columns,rows,meta):
        path=filedialog.asksaveasfilename(initialfile='异常检查.csv',defaultextension='.csv')
        if not path:return
        try:
            with open(path,'w',newline='',encoding='utf-8-sig') as f:
                writer=csv.writer(f);writer.writerow(['检查说明',json.dumps(meta,ensure_ascii=False)]);writer.writerow(columns);writer.writerows(rows)
        except Exception as e:self.error(e)
    def show_qc(self,record,draft=False):
        win=tk.Toplevel(self.root);win.title(record['plate']+' · 异常定位');win.geometry('1200x720')
        bar=ttk.Frame(win);bar.pack(fill='x',padx=12,pady=10)
        mode=tk.StringVar(value='当前阈值重查')
        options=['当前阈值重查']+(['保存时快照'] if record.get('qc_snapshot') and not draft else [])
        ttk.Combobox(bar,textvariable=mode,values=options,state='readonly',width=20).pack(side='left')
        info=ttk.Label(win,text='',wraplength=1150);info.pack(anchor='w',padx=12,pady=8)
        cols=['级别','检查项','孔位','时间 h','原因'];tr=table(win,cols,[85,180,180,100,600])
        tr.tag_configure('需核对',foreground='#a5372c')
        state={}
        def rebuild(*args):
            r=self.draft_record() if draft else record
            report=r['qc_snapshot'] if mode.get()=='保存时快照' else qc.inspect(r,self.qc_settings)
            state.update(record=r,report=report,rows=[],items={});tr.delete(*tr.get_children())
            for finding in report['findings']:
                row=[finding['level'],finding['title'],','.join(finding['wells']),','.join(map(str,finding['times'])),finding['reason']]
                iid=tr.insert('','end',values=row,tags=(finding['level'],));state['items'][iid]=finding;state['rows'].append(row)
            info.config(text=f"{mode.get()}：{report['alerts']} 条需核对，{report['unavailable']} 项无法评估。双击查看本孔、复孔和DMSO曲线。"+(' 修改孔位后点重查。' if draft else ' 保存时的结果与阈值在快照中保留。'))
        def selected():
            ids=tr.selection()
            return state['items'].get(ids[0]) if ids else None
        def curve():
            finding=selected()
            if not finding or not finding['wells']:return messagebox.showinfo('查看曲线','请选择含具体孔位的检查项',parent=win)
            self.detail(record=state['record'],focus=finding['wells'][0],highlight_times=finding['times'])
        def edit():
            finding=selected()
            if finding and finding['wells']:self.edit_well(finding['wells'][0])
        ttk.Button(bar,text='重查 / 刷新',command=rebuild).pack(side='left',padx=8)
        ttk.Button(bar,text='查看所选曲线',command=curve).pack(side='left')
        ttk.Button(bar,text='查看本次阈值',command=lambda:self.text_window('本次检查参数',json.dumps(state['report']['parameters'],ensure_ascii=False,indent=2)+'\n\n'+qc.HELP)).pack(side='left',padx=8)
        ttk.Button(bar,text='导出清单',command=lambda:self.export_qc(cols,state['rows'],dict(plate=state['record']['plate'],record_id=state['record'].get('id'),mode=mode.get(),parameters=state['report']['parameters'],qc_version=state['report']['version']))).pack(side='right')
        if draft:ttk.Button(bar,text='核对所选孔 / 手动排除',command=edit).pack(side='left')
        mode.trace_add('write',rebuild);tr.bind('<Double-1>',lambda e:curve());rebuild()
    def qc_options(self):
        win=tk.Toplevel(self.root);win.title('异常检查阈值');win.geometry('800x700');win.transient(self.root)
        ttk.Label(win,text='以下为探索性提示阈值，不是已验证标准。比例使用小数；不触发自动排除。',wraplength=750).pack(anchor='w',padx=15,pady=12)
        grid=ttk.Frame(win);grid.pack(fill='x',padx=15);fields={}
        for i,(key,label) in enumerate(qc.LABELS.items()):
            ttk.Label(grid,text=label).grid(row=i,column=0,sticky='w',pady=6)
            var=tk.StringVar(value='' if self.qc_settings[key] is None else str(self.qc_settings[key]));fields[key]=var
            ttk.Entry(grid,textvariable=var,width=14).grid(row=i,column=1,padx=20)
        def save():
            try:
                proposed=qc.config({k:v.get() for k,v in fields.items()})
                tmp=self.workspace/'qc_settings.tmp';tmp.write_text(json.dumps(proposed,ensure_ascii=False),encoding='utf-8');tmp.replace(self.workspace/'qc_settings.json')
                self.qc_settings=proposed;self.render_wells();self.refresh();win.destroy()
                self.status.set('异常检查参数已更新；历史保存时的快照保持不变')
            except Exception as e:self.error(e)
        bar=ttk.Frame(win);bar.pack(fill='x',padx=15,pady=15)
        ttk.Button(bar,text='规则与局限说明',command=lambda:self.text_window('异常检查规则',qc.HELP)).pack(side='left')
        ttk.Button(bar,text='保存并重查',command=save).pack(side='right')
    def refresh(self):
        self.history.delete(*self.history.get_children())
        for r in self.store.records(False):self.history.insert('','end',iid=str(r['id']),values=[r['id'],r['plate'],'有效' if r['active'] else '已移除 / 替换',r['created'],r['filename'],r.get('note','')])
        self.refresh_ranks()
        self.refresh_qc()
    def refresh_ranks(self):
        try:
            pars={k:float(v.get()) if v.get().strip() else None for k,v in self.params.items()}
            if any(v is not None and not math.isfinite(v) for v in pars.values()):raise ValueError('参数必须为有限数字')
            if pars['weight'] is None or not 0<=pars['weight']<=1:raise ValueError('权重必须在 0–1 之间')
            rec=self.store.records();rows=rankings(rec,**pars,qc_settings=self.qc_settings)
            rows=[r for r in rows if r['dV']>0 and (not self.mode.get().startswith('A') or r['H']<1)]
            sortkey={'综合分':'score','ΔV':'dV','ΔGR':'dGR','H GI':'HGI','H GR':'HGR','H Z':'HZ','H Robust Z':'HRZ','H 时间平均抑制':'HTI'}[self.sort.get()]
            lower=sortkey in ('HGR','HZ','HRZ')
            rows.sort(key=lambda r:(r[sortkey] is None,(r[sortkey] if lower else -r[sortkey]) if r[sortkey] is not None else 0,r['id']))
            self.ranktree.delete(*self.ranktree.get_children());self.rankrows=[]
            for i,r in enumerate(rows,1):
                values=[i]+[fmt(r[k]) for k in ['pair','id','name','score','EV','H','dV','dGR','HGI','HGR','HZ','HRZ','HTI','flags']]
                self.ranktree.insert('','end',values=values);self.rankrows.append(values)
            _,pending=pair_status(rec)
            self.rankinfo.config(text=f'{len(rec)} 块有效板 · {len(rows)} 条筛选结果 · 待配对：'+(', '.join(sorted(pending)) or '无')+'。门槛留空表示不限制；A 榜未设 EV 门槛时，不保证 EV 影响小。综合分是探索性排序规则，非统计显著性。')
            (self.workspace/'settings.json').write_text(json.dumps({k:v.get() for k,v in self.params.items()},ensure_ascii=False),encoding='utf-8')
            self.rank_snapshot={**{k:v.get() for k,v in self.params.items()},'mode':self.mode.get(),'sort':self.sort.get(),'records':[r['id'] for r in rec],'version':VERSION,'qc_settings':self.qc_settings}
            self.rank_data=rows
            if hasattr(self,'visual'):self.visual.update(rec,rows,self.rank_snapshot)
            if hasattr(self,'research'):self.research.update(rec,self.rank_snapshot)
        except Exception as e:self.error(e)
    def export_ranks(self):
        p=filedialog.asksaveasfilename(initialfile='药物排名.csv',defaultextension='.csv')
        if p:
            with open(p,'w',newline='',encoding='utf-8-sig') as f:
                writer=csv.writer(f);writer.writerow(['分析参数',json.dumps(self.rank_snapshot,ensure_ascii=False)]);writer.writerow(self.rankcols);writer.writerows(self.rankrows)
            self.status.set('已导出当前榜单和本次分析参数')

if __name__=='__main__':
    demo='--demo' in sys.argv
    if demo:ROOT=ROOT/'example'/'workspace'
    root=tk.Tk();app=App(root,demo=demo)
    if '--smoke' in sys.argv:root.after(700,app.close)
    root.mainloop()
