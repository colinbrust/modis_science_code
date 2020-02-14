import unittest
import ee
ee.Initialize()
from msc.models.GeeMod17 import MOD17
from msc.models.GeeMod16 import MOD16
from msc.utils.ExtractGeeData import ExtractGeeData
from msc.utils.Smap import get_smap


class MyTestCase(unittest.TestCase):

    def setUp(self) -> None:

        meteo = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')
        dayl = ee.ImageCollection('NASA/ORNL/DAYMET_V3').select('dayl')
        elev = ee.Image('USGS/NED')
        smap = get_smap(start='2003-01-01', end='2018-01-01')

        self.m17 = ExtractGeeData(template='../data/MOD17/template_example.csv',
                                  model=MOD17,
                                  out_dir='./tmp',
                                  meteorology=meteo)

        self.m16 = ExtractGeeData(template='../data/MOD17/template_example.csv',
                                  model=MOD16,
                                  out_dir='./tmp',
                                  meteorology=meteo,
                                  daylength=dayl,
                                  elev=elev,
                                  smap_sm=smap)

    def test_extraction(self):
        i = 0
        while i <= 5:
            try:
                test_df = self.m17._run_single_location(name='US-AR1', year=2010, lat=36.4267, lon=-99.42)
            except ee.ee_exception.EEException as e:
                print(e)
                i += 1
                continue
            break
        self.assertEqual(len(test_df), 365)

    def test_MOD16(self):
        i = 0
        while i <= 5:
            try:
                test_df = self.m16._run_single_location(name='US-AR1', year=2010, lat=36.4267, lon=-99.42)
            except ee.ee_exception.EEException as e:
                print(e)
                i += 1
                continue
            break
        self.assertIn('ET', test_df.columns)
        self.assertIn('name', test_df.columns)
        self.assertIn('system:time_start', test_df.columns)

    def test_MOD17(self):
        i = 0
        while i <= 5:
            try:
                test_df = self.m17._run_single_location(name='US-AR1', year=2010, lat=36.4267, lon=-99.42)
            except ee.ee_exception.EEException as e:
                print(e)
                i += 1
                continue
            break
        self.assertIn('GPP', test_df.columns)


if __name__ == '__main__':
    unittest.main()
