import unittest,copy,tempfile,json
from pathlib import Path
from engine import Store,parse_txt
from demo import synthetic_records,CURVES
import research,reviews

class ResearchTests(unittest.TestCase):
    def test_rule_categories_and_two_hour_onset(self):
        expected=['持续抑制','迟发抑制','抑制后恢复','持续抑制','持续抑制','未见持续抑制','短暂/波动','迟发抑制']
        self.assertEqual([research.classify(v)['label'] for v in CURVES],expected)
        self.assertEqual(research.classify(CURVES[1])['onset'],12)
        self.assertEqual(research.classify([1]*13)['TI'],0)
    def test_invalid_curve_is_not_imputed(self):
        values=[1]*13;values[5]=None
        self.assertEqual(research.classify(values)['label'],'无法分型')
        self.assertIsNone(values[5])
        for config in ({'consecutive':1},{'consecutive':2.5},{'inhibition':1,'recovery':.9},{'rebound':float('nan')}):
            with self.assertRaises(ValueError):research.temporal_config(config)
    def test_profiles_replicates_and_baseline_invariance(self):
        records=synthetic_records();before=copy.deepcopy(records);profiles=research.temporal_profiles(records)
        p=next(p for p in profiles if p['cell']=='R273H' and p['drug_id']=='DEMO-005')
        self.assertIn('不一致',p['agreement']);self.assertEqual(len(profiles),16);self.assertEqual(records,before)
        for t in research.TIMES:records[1]['data']['times'][str(t)]['A02']*=2
        again=research.temporal_profiles(records)
        p1=next(p for p in profiles if p['cell']=='R273H' and p['drug_id']=='DEMO-001')
        p2=next(p for p in again if p['cell']=='R273H' and p['drug_id']=='DEMO-001')
        self.assertEqual(p1['curve'],p2['curve'])
    def test_stability_reproducible_and_no_mutation(self):
        records=synthetic_records();before=copy.deepcopy(records)
        a=research.stability(records,top_k=3);b=research.stability(records,top_k=3)
        self.assertEqual(a,b);self.assertEqual(before,records)
        self.assertEqual(len(a['scenarios']),19)
        for r in a['summary']:
            self.assertEqual(r['eligible']+r['filtered']+r['unavailable'],r['total'])
            self.assertAlmostEqual(r['frequency'],r['top_hits']/r['total'])
        self.assertEqual(a['summary'][0]['top_hits'],19)
    def test_ties_no_arbitrary_topk_and_filtered_ranks(self):
        records=synthetic_records()
        # 令前两个药物的各组计数完全相同，使其分数形成真实并列。
        for r in records:
            for t in research.TIMES:
                for dst,src in [('A04','A02'),('A05','A03')]:r['data']['times'][str(t)][dst]=r['data']['times'][str(t)][src]
        out=research.stability(records,top_k=1)
        base=[d for d in out['details'] if d['scenario']==0 and d['drug_id'] in ('DEMO-001','DEMO-002')]
        self.assertTrue(all(d['rank']==1 and d['top'] for d in base))
        strict=research.stability(records,dict(ev_min=100),top_k=1)
        self.assertTrue(all(r['rank_min'] is None for r in strict['summary']))
    def test_unavailable_metric_and_missing_pair(self):
        records=synthetic_records()
        for r in records:
            for w in r['layout']:
                if w['role']=='DMSO对照':
                    for t in research.TIMES:r['data']['times'][str(t)][w['well']]=r['data']['times']['0'][w['well']]
        out=research.stability(records)
        self.assertTrue(any(r['unavailable']>0 for r in out['summary']))
        self.assertEqual(research.stability(records[:1])['summary'],[])
    def test_demo_raw_roundtrip(self):
        for r in synthetic_records():self.assertEqual(parse_txt(r['raw'])['times'],r['data']['times'])

class ReviewTests(unittest.TestCase):
    def test_append_history_require_evidence_and_source_changes(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as folder:
            s=Store(Path(folder)/'db.sqlite');records=synthetic_records()
            for r in records:s.save(r)
            sources=[dict(id=r['id'],plate=r['plate'],hash=r['hash']) for r in s.records()]
            d=dict(pair='01A',drug_id='DEMO-001',name='合成药物',status='待复筛',reviewer='测试员',notes='计划复筛',experiment_date='2026-10-30',independent_runs=0,source_records=sources)
            cid=reviews.append(s,d)
            self.assertEqual(reviews.source_state(reviews.latest(s)[0],s.records()),'源数据一致')
            invalid=dict(d,status='已复现')
            with self.assertRaises(ValueError):reviews.append(s,invalid,cid)
            update=dict(d,status='已复现',experiment_id='TEST-01',independent_runs=1,evidence='测试证据',revision_reason='追加测试结果')
            reviews.append(s,update,cid)
            events=reviews.history(s,cid);self.assertEqual(len(events),2);self.assertEqual(events[1]['status'],'待复筛')
            s.remove(s.records()[0]['id']);self.assertIn('变化',reviews.source_state(events[0],s.records()))
            run=research.stability(s.records());reviews.save_run(s,'stability',run);self.assertEqual(reviews.runs(s,'stability')[0]['result'],run)
            s.db.close();s=Store(Path(folder)/'db.sqlite');self.assertEqual(len(reviews.history(s,cid)),2);s.db.close()
    def test_evidence_hash_and_identity_guard(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as folder:
            s=Store(Path(folder)/'db.sqlite');evidence=Path(folder)/'evidence.txt';evidence.write_text('synthetic evidence')
            d=dict(pair='01A',drug_id='D1',name='合成',status='结果不确定',reviewer='测试员',notes='备注',experiment_date='2026-10-30',independent_runs=1,experiment_id='T1',evidence='测试描述',evidence_file=str(evidence))
            cid=reviews.append(s,d);event=reviews.latest(s)[0];self.assertEqual(len(event['evidence_sha256']),64)
            with self.assertRaises(ValueError):reviews.append(s,dict(d,drug_id='D2',revision_reason='改药物'),cid)
            self.assertEqual(len(reviews.history(s)),1);s.db.close()

if __name__=='__main__':unittest.main()
