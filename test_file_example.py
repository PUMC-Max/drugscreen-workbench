import copy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from engine import Store,TIMES,analyze,rankings
from file_example import map_rows,analysis_views,save_source,load_source
import research,reviews


def sample_rows():
    rows=[('plate_id','well','cell_line','treatment','replicate','treatment_class','include_in_plots','dmso_reference','notes')]
    named=['DMSO','RSL3','ML162','RSL3+Lip1','FIN56','iFSP1','erastin','Cells only']
    for i,row in enumerate('ABCDEFGH'):
        for col in range(1,13):
            within=(col-1)%6;cell='EV' if col<=6 else 'R273H'
            name=named[i] if within<2 else f'drug{i+1+(8 if within>=4 else 0)}'
            kind='Control/reference' if name in ('DMSO','Cells only') else 'Screen drug' if name.startswith('drug') else 'Named drug'
            rows.append((None,f'{row}{col}',cell,name,str(within%2+1),kind,False,name=='DMSO',None))
    return rows


def source():
    layout,_=map_rows(sample_rows());times={str(t):{} for t in TIMES}
    for w in layout:
        cell=w['cell'];n0=100 if cell=='EV' else 1000;end_fold=2 if cell=='EV' else 4
        for t in TIMES:
            effect=1 if w['role'] in ('DMSO对照','仅细胞') else 1-(.2 if cell=='EV' else .6)*t/24
            times[str(t)][w['well']]=n0*end_fold**(t/24)*effect
    return dict(plate='file-test',layout=layout,raw='source file unchanged',data=dict(times=times),provenance={})


class FileExampleTests(unittest.TestCase):
    def test_cell_specific_controls_and_all_treatments_without_mutation(self):
        s=source();before=copy.deepcopy(s);views=analysis_views(s);rows=rankings(views)
        self.assertEqual(len(rows),22);self.assertEqual(len(research.temporal_profiles(views)),44)
        for r in rows:
            self.assertAlmostEqual(r['EV'],.8);self.assertAlmostEqual(r['H'],.4);self.assertAlmostEqual(r['dV'],.4)
        self.assertEqual(s,before)
        for view in views:
            self.assertEqual(view['raw'],s['raw']);self.assertEqual(view['data'],s['data'])
            report=analyze(view)
            self.assertTrue(any(f['code'].startswith('positive_') and 'RSL3' in f['reason'] for f in report['qc']['findings']))
            self.assertIn('RSL3',{d['id'] for d in report['drugs']})
            self.assertTrue(all(not w['excluded'] for w in view['layout']))

    def test_plate_map_validation_and_plot_flag_is_not_exclusion(self):
        rows=sample_rows();layout,plate=map_rows(rows)
        self.assertIsNone(plate);self.assertEqual(sum(w['include_in_ranking'] for w in layout),88)
        self.assertTrue(all(not w['excluded'] for w in layout));self.assertEqual(next(w['dose'] for w in layout if w['name']=='ML162'),None)
        bad=copy.deepcopy(rows);bad[-1]=bad[1]
        with self.assertRaises(ValueError):map_rows(bad)
        bad=copy.deepcopy(rows);bad[1]=tuple(False if i==7 else v for i,v in enumerate(bad[1]))
        with self.assertRaises(ValueError):map_rows(bad)

    def test_source_snapshot_and_review_provenance(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as tmp:
            store=Store(Path(tmp)/'screen.sqlite')
            try:
                s=source();sid=save_source(store,s);self.assertEqual(sid,save_source(store,s))
                self.assertEqual(store.db.execute('SELECT COUNT(*) FROM file_example_sources').fetchone()[0],1)
                self.assertEqual(store.records(),[]);self.assertEqual(reviews.history(store),[])
                views=analysis_views(s,sid);event=dict(pair='file-test·',source_records=[dict(plate=r['plate'],id=r['id'],hash=r['hash']) for r in views])
                self.assertEqual(reviews.source_state(event,views),'源数据一致')
                s['data']['times']['24']['A03']*=2;changed=analysis_views(s,save_source(store,s))
                self.assertIn('变化',reviews.source_state(event,changed))
            finally:store.db.close()

    def test_only_plate_map_is_accessed(self):
        class Book:
            sheetnames=['Data','PlateMap'];closed=False
            def __getitem__(self,name):
                if name!='PlateMap':raise AssertionError('不允许读取 Excel Data')
                return type('Sheet',(),{'values':sample_rows()})()
            def close(self):self.closed=True
        book=Book();s=source()
        with patch('file_example.openpyxl.load_workbook',return_value=book),patch('file_example.read_txt',return_value=(s['raw'],s['data'])),patch.object(Path,'read_bytes',return_value=b'input'):
            loaded=load_source('layout.xlsx','counts.txt')
        self.assertEqual(loaded['data'],s['data']);self.assertTrue(book.closed)


if __name__=='__main__':unittest.main()
