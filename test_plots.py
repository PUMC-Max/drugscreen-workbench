import unittest,copy,tempfile,json,csv
from pathlib import Path
import plots
from matplotlib.backends.backend_agg import FigureCanvasAgg
from engine import analyze,rankings,Store
from test_qc import plate

def samples():
    e=plate();h=plate();h['plate']='01AH'
    for w in e['layout']:
        if w['role']=='DMSO对照':
            factor=1+int(w['well'][1:])*.001
            for r in (e,h):
                for t in plots.TIMES:r['data']['times'][str(t)][w['well']]*=factor
    for row in 'ABCDEFGH':
        for col in range(2,12):
            key=f'{row}{col:02}';end=.35+((ord(row)-65)*5+(col-2)//2)*.065
            for t in plots.TIMES:h['data']['times'][str(t)][key]=100*(1+(end-1)*t/24)
    h['data']['times']['12']['B02']=600
    e['raw']='synthetic EV';h['raw']='synthetic H'
    return [e,h]

class PlotTests(unittest.TestCase):
    def fig(self):
        f=plots.Figure(figsize=(11,6));FigureCanvasAgg(f);return f
    def test_heat_matches_engine_and_retains_flagged_value(self):
        r=samples()[1];a=analyze(r)
        for metric in ('V','GI','GR','Z','RZ'):
            rows=plots.heat_values(r,metric,24)
            for row in rows:
                if row['role']=='药物库':self.assertEqual(row['value'],a['metrics'][row['well']][metric])
        row=next(v for v in plots.heat_values(r,'Counts',12) if v['well']=='B02')
        self.assertTrue(row['qc_alert']);self.assertEqual(row['value'],600)
    def test_curves_reproduce_endpoint_and_exclusion(self):
        records=samples();r=records[1]
        for mode in ('V','相对0h','Counts'):
            data=plots.time_series(records,'01A','B-0',mode)
            self.assertEqual(len(data),52)
            if mode=='V':
                endpoint=[v['value'] for v in data if v['cell']=='R273H' and v['time']==24]
                self.assertAlmostEqual(sum(endpoint)/2,next(d for d in analyze(r)['drugs'] if d['id']=='B-0')['V'])
        next(w for w in r['layout'] if w['well']=='B02').update(excluded=True,reason='测试')
        result=plots.curves(self.fig(),records,'01A','B-0')
        raw=next(v for v in result['rows'] if v['cell']=='R273H' and v['well']=='B03' and v['time']==12)
        avg=next(v for v in result['rows'] if v['cell']=='R273H' and v['well']=='有效复孔均值' and v['time']==12)
        self.assertEqual(avg['value'],raw['value'])
        self.assertTrue(any(v['excluded'] for v in result['rows']))
    def test_empty_and_zero_and_unpaired(self):
        f=self.fig();self.assertEqual(plots.scatter(f,[])['rows'],[])
        self.assertEqual(plots.heatmap(f,None)['rows'],[])
        r=samples()[0];r['data']['times']['0']['B02']=0
        vals=plots.time_series([r],'01A','B-0','相对0h')
        self.assertTrue(all(v['value'] is None for v in vals if v['well']=='B02'))
        self.assertIn('尚缺数据',plots.curves(f,[r],'01A','B-0')['message'])
    def test_scatter_keeps_both_selectivity_directions(self):
        rows=rankings(samples());before=copy.deepcopy(rows)
        result=plots.scatter(self.fig(),rows)
        self.assertEqual(len(result['rows']),40)
        self.assertTrue(any(r['EV']<r['R273H'] for r in result['rows']))
        self.assertTrue(any(r['EV']>r['R273H'] for r in result['rows']))
        self.assertEqual(rows,before)
    def test_bar_sort_and_limit(self):
        rows=rankings(samples())
        for metric,key in plots.SORTS.items():
            result=plots.bars(self.fig(),rows,metric,10)
            values=[r['value'] for r in result['rows']]
            self.assertEqual(len(values),10)
            self.assertEqual(values,sorted(values,reverse=key not in ('HGR','HZ','HRZ')))
    def test_render_all_views_and_exports(self):
        records=samples();rows=rankings(records)
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as tmp:
            for name,fn in [('scatter',lambda f:plots.scatter(f,rows,'GR')),('bars',lambda f:plots.bars(f,rows,'ΔV')),
                            ('heat',lambda f:plots.heatmap(f,records[1],'V',12)),('curves',lambda f:plots.curves(f,records,'01A','B-0'))]:
                f=self.fig();data=fn(f);f.canvas.draw()
                for ext in ('png','svg'):
                    path=Path(tmp)/(name+'.'+ext);f.savefig(path,bbox_inches='tight');self.assertGreater(path.stat().st_size,1000)
                json.dumps(data['rows'],allow_nan=False)

if __name__=='__main__':unittest.main()
