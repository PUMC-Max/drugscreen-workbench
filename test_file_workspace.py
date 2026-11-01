import copy
import json
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch
import app
import file_example
import research
import reviews
from engine import analyze,rankings
from test_file_example import source


class FileWorkspaceTests(unittest.TestCase):
    def test_full_app_charts_edit_replace_remove_and_reopen(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as folder:
            out=Path(folder);s=source();s['filename']='example.txt';s['data']['header']='FILE EXAMPLE'
            with patch('file_example.load_source',return_value=s):workspace=file_example.prepare_workspace(destination=out)
            errors=[];root=app.tk.Tk();root.withdraw()
            with patch.object(app.App,'error',lambda self,e:errors.append(str(e))):
                a=file_example.open_workspace(root,workspace)
                try:
                    self.assertIs(type(a),app.App)
                    self.assertEqual({a.nb.tab(tab,'text') for tab in a.nb.tabs()},
                        {'新增 / 替换读板','实验记录','异常检查','两类药物榜单','结果可视化','方法与适用范围','时序分型','排名稳定性','复筛记录'})
                    self.assertEqual(a.nb.select(),str(a.vistab));self.assertEqual(len(a.history.get_children()),1)
                    v=a.visual;self.assertEqual(len(v.all),22);self.assertEqual(len(v.drug_map),22);self.assertEqual(len(v.result['rows']),22)
                    self.assertIn('1 块有效板',v.summary['text']);self.assertNotIn('待配对：file-test',v.summary['text'])
                    for kind in ('两组药效散点','选择性排名','96孔热图','单药时间曲线'):
                        v.kind.set(kind);v.render();v.canvas.draw()
                        self.assertTrue(v.result['rows'],kind);self.assertNotIn('暂未生成',v.result['message'])
                    v.kind.set('96孔热图');v.heatmetric.set('V');v.time.set('24');v.render()
                    heat={r['well']:r for r in v.result['rows']};self.assertEqual(len(heat),96)
                    self.assertAlmostEqual(heat['A03']['value'],.8);self.assertAlmostEqual(heat['A09']['value'],.4)
                    v.heatmetric.set('GR');v.render();metrics=analyze(a.store.records()[0])['metrics']
                    self.assertAlmostEqual(v.result['heat']['A09']['value'],metrics['A09']['GR'])
                    v.focus('file-test·','RSL3');self.assertEqual(len(v.result['rows']),78)
                    self.assertEqual({r['cell'] for r in v.result['rows']},{'EV','R273H'})
                    with patch('visualization.filedialog.asksaveasfilename',return_value=str(out/'chart.png')):v.export()
                    self.assertTrue((out/'chart.png').exists());meta=json.loads((out/'chart.meta.json').read_text(encoding='utf-8'))
                    self.assertEqual(len(meta['records']),1)
                    # 缺失浓度保持空值；可以用与正式版相同的板位编辑器保存。
                    a.edit_well('C01');win=next(w for w in root.winfo_children() if isinstance(w,app.tk.Toplevel))
                    next(w for w in win.winfo_children() if isinstance(w,app.ttk.Button) and w['text']=='保存本孔修改').invoke()
                    self.assertIsNone(next(w['dose'] for w in a.layout if w['well']=='C01'))
                    panel=a.research;self.assertEqual(len(panel.profiles),44);panel.start_stability();deadline=time.monotonic()+15
                    while panel.busy and time.monotonic()<deadline:root.update();time.sleep(.02)
                    self.assertFalse(panel.busy);self.assertEqual(len(panel.result['summary']),22)
                    event=dict(pair='file-test·',drug_id='drug1',name='drug1',status='待复筛',reviewer='test',notes='plan',experiment_date='2026-10-30',independent_runs=0,
                               source_records=[dict(id=r['id'],plate=r['plate'],hash=r['hash']) for r in panel.records])
                    case=reviews.append(a.store,event);self.assertEqual(reviews.source_state(event,panel.records),'源数据一致')
                    changed=copy.deepcopy(a.store.records()[0]);changed['raw']+='\nreplacement';changed['data']['times']['24']['A03']*=.5
                    a.store.save(changed,replace=True);a.refresh()
                    self.assertEqual(len(a.store.records()),1);self.assertEqual(len(a.history.get_children()),2)
                    self.assertIn('过期',panel.stability_info['text']);self.assertIn('变化',reviews.source_state(event,panel.records))
                    a.history.selection_set(str(a.store.records()[0]['id']))
                    with patch('app.messagebox.askyesno',return_value=True):a.remove()
                    self.assertEqual(a.store.records(),[]);self.assertEqual(a.rank_data,[]);self.assertEqual(panel.profiles,[])
                    v.kind.set('两组药效散点');v.render();self.assertEqual(v.result['rows'],[])
                    self.assertEqual(len(reviews.history(a.store,case)),1)
                    with patch('file_example.load_source') as load:file_example.prepare_workspace(destination=out);load.assert_not_called()
                    self.assertEqual(a.store.records(),[]);root.update();self.assertFalse(errors,errors)
                finally:a.close()

    def test_physical_record_matches_derived_views(self):
        s=source();views=file_example.analysis_views(s)
        self.assertEqual(rankings([s]),rankings(views))
        self.assertEqual(len(research.temporal_profiles([s])),44)
        report=analyze(s)
        self.assertEqual(len(report['drugs']),44);self.assertEqual(len(report['metrics']),96)
        self.assertEqual(report['control_folds'],{'EV':2.,'R273H':4.})
        for f in report['qc']['findings']:
            self.assertIn(f['cell'],('EV','R273H'))


if __name__=='__main__':unittest.main()
