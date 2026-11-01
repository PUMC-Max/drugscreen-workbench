"""物理板与细胞分析分组。物理记录保留96孔，分组只用于计算。"""
import copy
import hashlib
import json


def cell_name(well):
    value=well.get('cell')
    return 'R273H' if value=='H' else value


def mixed_cells(record):
    if record.get('view_only'):return False
    return {'EV','R273H'}.issubset({cell_name(w) for w in record['layout']})


def cell_views(record):
    if not mixed_cells(record):return [record]
    result=[]
    for cell,suffix in [('EV','E'),('R273H','H')]:
        layout=copy.deepcopy(record['layout'])
        for w in layout:
            w['source_role']=w['role']
            if cell_name(w)!=cell:w.update(role='未使用',include_in_ranking=False)
        view=dict(record,layout=layout,plate=record['plate']+'·'+suffix,physical_plate=record['plate'],
            physical_record_id=record.get('id'),analysis_cell=cell,id=f"{record.get('id')}:{cell}",view_only=True)
        view['hash']=hashlib.sha256(json.dumps(dict(layout=layout,data=record['data']),sort_keys=True,ensure_ascii=False).encode()).hexdigest()
        result.append(view)
    return result


def analysis_records(records):
    return [view for record in records for view in cell_views(record)]


def pair_status(records):
    plates={r['plate'] for r in analysis_records(records)}
    paired=sum(p.endswith('E') and p[:-1]+'H' in plates for p in plates)
    pending=sorted(p for p in plates if p[:-1]+('H' if p.endswith('E') else 'E') not in plates)
    return paired,pending
