import unittest, tempfile, copy
from pathlib import Path
from engine import *

def fixture(plate='01AE',drugfold=2):
    layout=[];data={str(t):{} for t in TIMES}
    for row in 'ABCDEFGH':
        for col in range(1,13):
            key=f'{row}{col:02}'
            role='DMSO对照' if col==1 else '药物库' if row=='A' and col in (2,3) else '未使用'
            layout.append(dict(well=key,role=role,drug_id='DRUG' if role=='药物库' else '',name='药物',dose=5 if role=='药物库' else 0,excluded=False))
            for t in TIMES:data[str(t)][key]=100*(1+(drugfold-1)*t/24) if role=='药物库' else 100*(1+t/24)
    lines=['Metric: Phase Object Count (Per Image)']
    for t in TIMES:
        lines.extend([f'Time Stamp:\ttest\tElapsed:\t{t}\thours','\t'+'\t'.join(str(x) for x in range(1,13))])
        for row in 'ABCDEFGH':lines.append(row+'\t'+'\t'.join(str(data[str(t)][f'{row}{col:02}']) for col in range(1,13)))
    return dict(plate=plate,layout=layout,data={'times':data},raw='\n'.join(lines),filename='test.txt')

class Tests(unittest.TestCase):
    def test_no_effect_and_arrest(self):
        d=analyze(fixture())['drugs'][0]
        self.assertAlmostEqual(d['V'],1);self.assertAlmostEqual(d['GI'],0);self.assertAlmostEqual(d['GR'],1);self.assertAlmostEqual(d['TI'],0)
        d=analyze(fixture(drugfold=1))['drugs'][0]
        self.assertAlmostEqual(d['GR'],0);self.assertAlmostEqual(d['GI'],1);self.assertAlmostEqual(d['V'],.5)
    def test_death_limit_and_invalid_control(self):
        d=analyze(fixture(drugfold=0))['drugs'][0];self.assertEqual(d['GR'],-1);self.assertEqual(d['GI'],2)
        r=fixture()
        for t in TIMES:
            for row in 'ABCDEFGH':r['data']['times'][str(t)][row+'01']=100
        self.assertIsNone(analyze(r)['drugs'][0]['GR'])
    def test_zero_baseline_invalidates_aggregate(self):
        r=fixture();r['data']['times']['0']['A02']=0
        self.assertIsNone(analyze(r)['drugs'][0]['GR'])
    def test_parser_rejects_bad_data(self):
        r=fixture();self.assertEqual(len(parse_txt(r['raw'])['times']),13)
        for bad in [r['raw']+'\n'+r['raw'],r['raw'].replace('100.0','NaN',1),r['raw'].replace('A\t100.0','A\t',1),r['raw'].replace('Elapsed:\t24','Elapsed:\t26')]:
            with self.assertRaises(ValueError):parse_txt(bad)
    def test_pairing_filters_and_direction(self):
        e=fixture();h=fixture('01AH',1)
        self.assertEqual(rankings([e]),[])
        result=rankings([e,h])[0];self.assertEqual(result['dV'],.5);self.assertEqual(result['dGR'],1)
        self.assertEqual(rankings([e,h],ev_min=1.1),[])
        self.assertEqual(rankings([e,fixture('02AH',1)]),[])
    def test_replace_remove_reopen_and_backup(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as tmp:
            p=Path(tmp)/'db.sqlite';s=Store(p);e=fixture();h=fixture('01AH',1)
            s.save(e);s.save(h);self.assertEqual(len(rankings(s.records())),1)
            with self.assertRaises(ValueError):s.save(e)
            s.save(fixture(drugfold=3),replace=True);self.assertEqual(len(s.records(False)),3)
            rid=next(r['id'] for r in s.records() if r['plate']=='01AH');s.remove(rid)
            self.assertEqual(rankings(s.records()),[]);s.db.close()
            reopened=Store(p);self.assertEqual(len(reopened.records()),1);self.assertEqual(len(reopened.records(False)),3)
            self.assertTrue(list((Path(tmp)/'backups').glob('*.sqlite')));reopened.db.close()
    def test_duplicate_file_across_plates(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as tmp:
            s=Store(Path(tmp)/'db');s.save(fixture());r=fixture('01AH')
            with self.assertRaises(ValueError):s.save(r)
            s.db.close()

if __name__=='__main__':unittest.main()
