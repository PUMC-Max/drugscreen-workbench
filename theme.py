"""共享桌面外观与侧边导航，所有工作区使用同一界面。"""
import tkinter as tk
from tkinter import ttk

def apply(root):
    root.configure(bg='#f3f6fa')
    s=ttk.Style(root);s.theme_use('clam')
    s.configure('.',font=('Microsoft YaHei UI',10),background='#f3f6fa',foreground='#24364b')
    s.configure('TFrame',background='#f3f6fa')
    s.configure('TLabel',background='#f3f6fa',foreground='#24364b')
    s.configure('TButton',background='#ffffff',foreground='#314b64',bordercolor='#dde5ee',lightcolor='#ffffff',darkcolor='#ffffff',padding=(13,9),relief='flat')
    s.map('TButton',background=[('active','#e8f1fa'),('pressed','#d7e8f7')])
    s.configure('TCombobox',fieldbackground='white',background='white',arrowcolor='#657b91',padding=6)
    s.map('TCombobox',fieldbackground=[('readonly','white')],selectbackground=[('readonly','white')],selectforeground=[('readonly','#24364b')])
    s.configure('TEntry',fieldbackground='white',padding=7,bordercolor='#dde5ee')
    s.configure('Treeview',background='white',fieldbackground='white',rowheight=34,borderwidth=0)
    s.configure('Treeview.Heading',background='#eaf0f7',foreground='#52677d',padding=9,relief='flat',font=('Microsoft YaHei UI',10,'bold'))
    s.map('Treeview',background=[('selected','#dcecf9')],foreground=[('selected','#17518a')])
    s.configure('Title.TLabel',font=('Microsoft YaHei UI',23,'bold'),foreground='#172c40')
    s.configure('Workspace.TNotebook',background='#f3f6fa',borderwidth=0,tabmargins=0)
    s.layout('Workspace.TNotebook.Tab',[])

def navigation(parent,notebook):
    tk.Label(parent,text='SCREEN / LAB',bg='#172c40',fg='#8da5bd',font=('Segoe UI',10,'bold'),anchor='w').pack(fill='x',padx=20,pady=(25,18))
    buttons={}
    for tab in notebook.tabs():
        button=tk.Button(parent,text=notebook.tab(tab,'text'),anchor='w',relief='flat',bd=0,bg='#172c40',fg='#b9c9d9',activebackground='#284966',activeforeground='white',
            font=('Microsoft YaHei UI',10),padx=18,pady=12,command=lambda tab=tab:notebook.select(tab),cursor='hand2')
        button.pack(fill='x',padx=10,pady=3);buttons[tab]=button
    def update(*_):
        for tab,button in buttons.items():button.configure(bg='#2d5679' if notebook.select()==tab else '#172c40',fg='white' if notebook.select()==tab else '#b9c9d9')
    notebook.bind('<<NotebookTabChanged>>',update,add='+');update()
    tk.Label(parent,text='本地实验工作区\n板位 · 数据 · 分析 · 复筛',justify='left',bg='#172c40',fg='#7e96ad',font=('Microsoft YaHei UI',9)).pack(side='bottom',anchor='w',padx=20,pady=24)
