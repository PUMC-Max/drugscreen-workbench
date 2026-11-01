import openpyxl
from pathlib import Path
from engine import validate_layout

def load_excel(path,default_plate=None):
    wb=openpyxl.load_workbook(path,read_only=True,data_only=True)
    try:
        if '实验板逐孔信息' not in wb.sheetnames:
            if 'PlateMap' not in wb.sheetnames:raise ValueError('需要“实验板逐孔信息”或 PlateMap 工作表')
            from file_example import map_rows
            layout,plate=map_rows(wb['PlateMap'].values)
            return {plate or default_plate or Path(path).stem:layout}
        rows=iter(wb['实验板逐孔信息'].values);header=next(rows)
        required=['样本板号','细胞','样本孔位','孔类型','MCE编号','药物名称']
        if any(k not in header for k in required):raise ValueError('Excel 缺少所需列：'+', '.join(required))
        output={}
        for row in rows:
            d=dict(zip(header,row));well=d.get('样本孔位')
            if not well:continue
            plate=str(d['样本板号']).strip();cell=d['细胞']
            if not ((plate.endswith('E') and cell=='EV') or (plate.endswith('H') and cell in ('H','R273H'))):
                raise ValueError(f'{plate}：细胞类型与板号后缀不符')
            role=d['孔类型']
            output.setdefault(plate,[]).append(dict(well=str(well).strip(),role=role,
                drug_id=str(d.get('MCE编号') or ''),name=str(d.get('药物名称') or ''),
                dose=5 if role=='药物库' else 1 if role=='阳性药物对照' else 0,
                excluded=False,reason='',source_plate=d.get('原库板号'),source_well=d.get('原库孔位'),
                replicate=d.get('复孔编号'),cell=cell))
        if not output:raise ValueError('未读取到板位')
        for layout in output.values():validate_layout(layout)
        return output
    finally:wb.close()
