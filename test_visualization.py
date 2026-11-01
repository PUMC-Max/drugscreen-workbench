import unittest,tempfile,json,csv
from pathlib import Path
from unittest.mock import patch
import app
from test_plots import samples

class VisualizationUITests(unittest.TestCase):
    def test_refresh_drilldown_export_and_remove(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as folder:
            base=Path(folder);records=samples()
            (base/'templates.json').write_text(json.dumps({r['plate']:r['layout'] for r in records}),encoding='utf-8')
            errors=[];root=app.tk.Tk();root.withdraw()
            with patch.object(app,'ROOT',base),patch.object(app.App,'error',lambda self,e:errors.append(str(e))):
                a=app.App(root)
                try:
                    for r in records:a.store.save(r)
                    a.refresh();v=a.visual;self.assertEqual(len(v.result['rows']),40)
                    v.scope.set('当前筛选榜单');v.render();self.assertEqual(len(v.result['rows']),len(a.rank_data))
                    a.params['ev_min'].set('99');a.refresh_ranks();self.assertEqual(v.result['rows'],[])
                    v.scope.set('全部已配对药物');v.render();self.assertEqual(len(v.result['rows']),40)
                    path=base/'figure.png'
                    with patch('visualization.filedialog.asksaveasfilename',return_value=str(path)):v.export()
                    self.assertTrue(path.exists())
                    meta=json.loads(path.with_suffix('.meta.json').read_text(encoding='utf-8'))
                    self.assertEqual(meta['scope'],'全部已配对药物');self.assertEqual(len(meta['records']),2)
                    with path.with_suffix('.data.csv').open(encoding='utf-8-sig') as f:self.assertEqual(len(list(csv.DictReader(f))),40)
                    v.focus('01A','B-0');self.assertEqual(v.kind.get(),'单药时间曲线');self.assertEqual(len(v.result['rows']),78)
                    v.kind.set('96孔热图');v.plate.set('01AH');v.render();self.assertEqual(len(v.result['rows']),96)
                    hid=next(r['id'] for r in a.store.records() if r['plate']=='01AH');a.store.remove(hid);a.refresh()
                    self.assertNotIn('01AH',v.platebox['values'])
                    v.kind.set('两组药效散点');v.render();self.assertEqual(v.result['rows'],[])
                    root.update();self.assertFalse(errors,errors)
                finally:a.store.db.close();root.destroy()

if __name__=='__main__':unittest.main()
