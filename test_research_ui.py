import unittest,tempfile,json,time
from pathlib import Path
from unittest.mock import patch
import app,reviews
from demo import synthetic_records

class ResearchUITests(unittest.TestCase):
    def test_background_run_snapshot_and_stale_notice(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as folder:
            base=Path(folder);records=synthetic_records()
            (base/'templates.json').write_text(json.dumps({r['plate']:r['layout'] for r in records}),encoding='utf-8')
            root=app.tk.Tk();root.withdraw();errors=[]
            with patch.object(app,'ROOT',base),patch.object(app.App,'error',lambda self,e:errors.append(str(e))):
                a=app.App(root,demo=True)
                try:
                    for r in records:a.store.save(r)
                    a.refresh();panel=a.research;self.assertEqual(len(panel.profiles),16)
                    panel.start_stability();deadline=time.monotonic()+15
                    while panel.busy and time.monotonic()<deadline:root.update();time.sleep(.02)
                    self.assertFalse(panel.busy);self.assertEqual(len(panel.result['summary']),8)
                    self.assertIn('一致',panel.stability_info['text']);self.assertEqual(len(reviews.runs(a.store,'stability')),1)
                    panel.temporal_tree.selection_set('0');panel.profile_plot();panel.temporal_options();panel.review_from_profile()
                    # 改变当前条件使旧快照过期，旧结果仍完整保存。
                    a.params['ev_min'].set('.85');a.refresh_ranks();self.assertIn('过期',panel.stability_info['text'])
                    root.update();self.assertFalse(errors,errors)
                finally:a.store.db.close();root.destroy()

if __name__=='__main__':unittest.main()
