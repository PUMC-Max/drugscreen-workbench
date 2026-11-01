import unittest,copy
import overview_plots as overview
import plots
from matplotlib.backends.backend_agg import FigureCanvasAgg

class OverviewTests(unittest.TestCase):
    def test_units_signs_and_missing(self):
        r=dict(EVGI=-.1,HGI=.2,EV=1.1,H=.8,EVTI=-.1,HTI=.2,EVGR=1.1,HGR=.8)
        for m in overview.METRICS[:3]:
            ev,h=overview.effects(r,m);self.assertAlmostEqual(ev,-10);self.assertAlmostEqual(h,20)
        self.assertEqual(overview.effects(dict(r,HGI=None),overview.METRICS[0]),None)
        ev,h=overview.effects(r,overview.METRICS[3]);self.assertAlmostEqual(h-ev,.3)
    def test_full_population_and_paging_keep_negative_values(self):
        rows=[dict(pair='01A',id=str(i),name='drug'+str(i),EVGI=.5,HGI=i/100,flags='需核对') for i in range(70)]
        before=copy.deepcopy(rows);fig=plots.Figure(figsize=(12,8));FigureCanvasAgg(fig)
        first=overview.waterfall(fig,rows);second=overview.waterfall(fig,rows,offset=30);last=overview.waterfall(fig,rows,offset=60)
        self.assertEqual(len(first['rows'])+len(second['rows'])+len(last['rows']),70)
        self.assertGreater(first['rows'][0]['value'],0);self.assertLess(last['rows'][-1]['value'],0)
        self.assertEqual(last['rows'][-1]['rank'],70);self.assertEqual(rows,before)
        result=overview.scatter(fig,rows);fig.canvas.draw();self.assertEqual(len(result['rows']),70)
        self.assertEqual(result['rows'][0]['difference'],-50)

if __name__=='__main__':unittest.main()
