"""复筛为人工提交的独立记录，采用仅追加的版本历史。"""
import json, hashlib
from datetime import date, datetime, timezone, timedelta

STATUSES=['待复筛','复筛中','已复现','未复现','结果不确定','暂不推进','记录作废']

def initialize(store):
    store.db.execute('CREATE TABLE IF NOT EXISTS rescreen_events (id INTEGER PRIMARY KEY, case_id INTEGER NOT NULL, created TEXT NOT NULL, payload TEXT NOT NULL)')
    store.db.execute('CREATE TABLE IF NOT EXISTS analysis_runs (id INTEGER PRIMARY KEY, created TEXT NOT NULL, kind TEXT NOT NULL, payload TEXT NOT NULL)')
    store.db.commit()

def now():return datetime.now(timezone(timedelta(hours=8))).isoformat(timespec='seconds')

def append(store,data,case_id=None):
    initialize(store);d=dict(data)
    for key in ('pair','drug_id','name','status','reviewer','notes'):
        if not str(d.get(key,'')).strip():raise ValueError('请填写 '+{'pair':'配对编号','drug_id':'药物编号','name':'药物名称','status':'状态','reviewer':'记录人','notes':'记录内容'}[key])
    if d['status'] not in STATUSES:raise ValueError('未知复筛状态')
    date.fromisoformat(d['experiment_date'])
    value=d.get('independent_runs',0)
    if isinstance(value,bool) or str(value).strip()!=str(int(value)) or int(value)<0:raise ValueError('独立实验次数需为非负整数')
    d['independent_runs']=int(value)
    if d['status'] in ('已复现','未复现','结果不确定') and (d['independent_runs']<1 or not d.get('experiment_id') or not d.get('evidence')):
        raise ValueError('提交复筛结果需填写实验编号、至少一次独立实验和证据说明；技术复孔不计为独立实验')
    if case_id is not None:
        old=next((e for e in latest(store) if e['case_id']==case_id),None)
        if not old:raise ValueError('复筛记录不存在')
        if (old['pair'],old['drug_id'])!=(d['pair'],d['drug_id']):raise ValueError('版本不能变更药物身份；请新建记录')
        if not str(d.get('revision_reason','')).strip():raise ValueError('更新记录需填写版本变更原因')
    file=d.get('evidence_file','')
    if file:
        sha=hashlib.sha256()
        try:
            with open(file,'rb') as f:
                for block in iter(lambda:f.read(1024*1024),b''):sha.update(block)
        except OSError:raise ValueError('证据文件不可读取；请选择存在的文件或仅填写证据说明')
        d['evidence_sha256']=sha.hexdigest()
    else:d['evidence_sha256']=None
    d['label_origin']='人工复筛记录（未经系统自动验证）'
    store.backup()
    with store.db:
        if case_id is None:case_id=store.db.execute('SELECT COALESCE(MAX(case_id),0)+1 FROM rescreen_events').fetchone()[0]
        store.db.execute('INSERT INTO rescreen_events(case_id,created,payload) VALUES(?,?,?)',(case_id,now(),json.dumps(d,ensure_ascii=False,allow_nan=False)))
    return case_id

def history(store,case_id=None):
    initialize(store)
    sql='SELECT id,case_id,created,payload FROM rescreen_events'
    params=()
    if case_id is not None:sql+=' WHERE case_id=?';params=(case_id,)
    return [dict(json.loads(p),event_id=i,case_id=c,created=t) for i,c,t,p in store.db.execute(sql+' ORDER BY id DESC',params)]

def latest(store):
    out={}
    for e in history(store):out.setdefault(e['case_id'],e)
    return list(out.values())

def source_state(event,records):
    old={(r['plate'],r.get('id'),r.get('hash')) for r in event.get('source_records',[])}
    current={(r['plate'],r.get('id'),r.get('hash')) for r in records if r['plate'][:-1]==event['pair']}
    return '源数据一致' if old and old==current else '源数据已变化 / 缺失，需复核'

def save_run(store,kind,result):
    initialize(store);store.backup()
    with store.db:
        cur=store.db.execute('INSERT INTO analysis_runs(created,kind,payload) VALUES(?,?,?)',(now(),kind,json.dumps(result,ensure_ascii=False,allow_nan=False)))
    return cur.lastrowid

def runs(store,kind):
    initialize(store)
    return [dict(id=i,created=t,result=json.loads(p)) for i,t,p in store.db.execute('SELECT id,created,payload FROM analysis_runs WHERE kind=? ORDER BY id DESC',(kind,))]
