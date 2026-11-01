"""Portable entry point. Default: isolated synthetic demo."""
import argparse,json
from pathlib import Path
import tkinter as tk
from app import App
from demo import generate
from templates import load_excel

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--layout',type=Path,help='Create/update the formal layout from your own Excel file')
    parser.add_argument('--formal',action='store_true',help='Open the formal local workspace')
    parser.add_argument('--smoke',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args();here=Path(__file__).resolve().parent
    if args.layout:
        (here/'templates.json').write_text(json.dumps(load_excel(args.layout),ensure_ascii=False),encoding='utf-8')
    formal=args.formal or bool(args.layout)
    if formal:
        workspace=here
        if not (workspace/'templates.json').exists():parser.error('Provide --layout YOUR_LAYOUT.xlsx for the first formal launch.')
    else:
        workspace=here/'example'/'workspace'
        if not (workspace/'templates.json').exists():generate()
    root=tk.Tk();app=App(root,demo=not formal,workspace=workspace)
    if args.smoke:root.after(700,app.close)
    root.mainloop()

if __name__=='__main__':main()
