import copy
import tempfile
import unittest
from pathlib import Path
import qc
from engine import analyze, rankings, Store, TIMES
from test_engine import fixture

def plate():
    r=fixture()
    for w in r['layout']:
        key=w['well'];row=key[0];col=int(key[1:])
        w.update(role='药物库',drug_id=f'{row}-{(col-2)//2}',name='药物',dose=5)
        if col in (1,12):
            w['role']='DMSO对照' if row in 'ABC' else '阳性药物对照' if row in 'DEF' else '仅细胞'
            w['name']='RSL3' if row in 'DEF' else '';w['dose']=1 if row in 'DEF' else 0
        for t in TIMES:r['data']['times'][str(t)][key]=100*(1+t/24) if w['role']!='阳性药物对照' else 100
    return r

def codes(r,settings=None):return {f['code'] for f in qc.inspect(r,settings)['findings']}

class QualityTests(unittest.TestCase):
    def test_clean_plate_and_no_mutation(self):
        r=plate();before=copy.deepcopy(r);report=qc.inspect(r)
        self.assertEqual(report['alerts'],0);self.assertEqual(report['unavailable'],0);self.assertEqual(before,r)
    def test_spike_localization(self):
        r=plate();r['data']['times']['12']['B02']=600
        spike=next(f for f in qc.inspect(r)['findings'] if f['code']=='spike')
        self.assertEqual(spike['wells'],['B02']);self.assertEqual(spike['times'],[12])
    def test_same_endpoint_different_curve(self):
        r=plate()
        for t in (6,8,10,12,14,16,18):r['data']['times'][str(t)]['B02']*=2
        self.assertIn('replicate_curve',codes(r));self.assertNotIn('replicate_end',codes(r))
    def test_initial_outlier_and_zero(self):
        r=plate();r['data']['times']['0']['B02']=400;r['data']['times']['0']['C02']=0
        self.assertTrue({'baseline','zero_start','replicate_unavailable'}<=codes(r))
    def test_control_leave_one_out(self):
        r=plate();r['data']['times']['12']['A01']=800
        hits=[f for f in qc.inspect(r)['findings'] if f['code']=='control_outlier']
        self.assertEqual(len(hits),1);self.assertEqual(hits[0]['wells'],['A01']);self.assertEqual(hits[0]['times'],[12])
        self.assertIn('control_cv',codes(r))
    def test_excluded_control_does_not_affect_reference(self):
        r=plate();r['data']['times']['12']['A01']=800
        next(w for w in r['layout'] if w['well']=='A01').update(excluded=True,reason='拍摄失焦')
        self.assertNotIn('control_outlier',codes(r));self.assertNotIn('control_cv',codes(r));self.assertIn('excluded',codes(r))
    def test_positive_response_and_disabled_gate(self):
        r=plate()
        for w in r['layout']:
            if w['role']=='阳性药物对照':r['data']['times']['24'][w['well']]=200
        self.assertIn('positive_weak',codes(r));self.assertNotIn('positive_weak',codes(r,{'positive_max':None}))
        for w in r['layout']:
            if w['role']=='阳性药物对照':w['excluded']=True
        self.assertIn('positive_missing',codes(r))
    def test_spatial_baseline_only(self):
        r=plate()
        for c in range(2,12):r['data']['times']['24'][f'B{c:02}']=400
        self.assertNotIn('spatial',codes(r))
        for c in range(2,12):r['data']['times']['0'][f'B{c:02}']=200
        hit=next(f for f in qc.inspect(r)['findings'] if f['code']=='spatial')
        self.assertIn('B行',hit['reason']);self.assertEqual(hit['times'],[0])
    def test_thresholds_and_invalid_parameters(self):
        r=plate();r['data']['times']['12']['B02']=600
        self.assertNotIn('spike',codes(r,{'spike_rel':10}))
        for value in ('x',-1,0,float('nan'),float('inf')):
            with self.assertRaises(ValueError):qc.config({'pair_rms':value})
    def test_real_strong_smooth_effect_not_removed(self):
        r=plate()
        for k in ('B02','B03'):
            for t in TIMES:r['data']['times'][str(t)][k]=100*.25**(t/24)
        self.assertEqual(qc.inspect(r)['alerts'],0)
        result=analyze(r);drug=next(d for d in result['drugs'] if d['wells']=='B02,B03')
        self.assertLess(drug['GR'],0);self.assertEqual(drug['n'],2)
    def test_ranking_warnings_do_not_change_effects(self):
        e=plate();h=plate();h['plate']='01AH'
        h['data']['times']['12']['B02']=600
        loose=rankings([e,h],qc_settings={'spike_rel':10,'pair_rms':10,'jump_rel':10})
        strict=rankings([e,h]);a=next(r for r in loose if r['id']=='B-0');b=next(r for r in strict if r['id']=='B-0')
        self.assertEqual(a['dV'],b['dV']);self.assertEqual(a['dGR'],b['dGR']);self.assertIn('孤立尖峰',b['flags'])
    def test_all_computable_metrics_survive_strict_quality_flags(self):
        r=plate()
        # 给对照微小变异，使普通和稳健Z均可计算。
        for i,k in enumerate(('A01','A12','B01','B12','C01','C12')):
            for t in TIMES:r['data']['times'][str(t)][k]*=1+i*.01
        r['data']['times']['12']['B02']=600
        strict={k:.0001 for k in qc.DEFAULTS};strict['min_growth']=10
        before=analyze(r);after=analyze(r,strict)
        self.assertGreater(after['qc']['alerts'],before['qc']['alerts'])
        self.assertEqual(before['metrics'],after['metrics'])
        for a,b in zip(before['drugs'],after['drugs']):
            for key in ('V','GI','GR','Z','RZ','TI'):
                self.assertIsNotNone(a[key]);self.assertEqual(a[key],b[key])
        self.assertEqual(sum(w.get('excluded',False) for w in r['layout']),0)
    def test_snapshot_and_history_recheck(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as tmp:
            s=Store(Path(tmp)/'test.sqlite');r=plate();r['data']['times']['12']['B02']=600;r['qc_settings']={'spike_rel':10}
            s.save(r);saved=s.records()[0];self.assertEqual(saved['qc_snapshot']['parameters']['spike_rel'],10)
            qc.inspect(saved,{'spike_rel':.1})
            self.assertEqual(s.records()[0]['qc_snapshot'],saved['qc_snapshot']);s.db.close()

if __name__=='__main__':unittest.main()
