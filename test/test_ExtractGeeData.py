import unittest
import ee
from msc.models.GeeMod17 import MOD17
from msc.utils.ExtractGeeData import ExtractGeeData

ee.Initialize()


class MyTestCase(unittest.TestCase):

    def setUp(self) -> None:

        meteo = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')
        self.extractor = ExtractGeeData(template='../data/MOD17/template_example.csv',
                                        model=MOD17,
                                        out_dir='',
                                        meteorology=meteo)

    def test_model_extraction(self):

        test_df = self.extractor._run_single_location(name='US-AR1', year=2010, lat=36.4267, lon=-99.42)
        self.assertEqual(len(test_df), 365)


if __name__ == '__main__':
    unittest.main()
