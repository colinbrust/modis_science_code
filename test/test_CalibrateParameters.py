import unittest
import pandas as pd
from msc.models.LocalMod17 import MOD17
from msc.utils import DatasetUtils as du
from msc.calval.CalibrateParameters import CalibrateParameters
import glob
import os
import numpy as np


class MyTestCase(unittest.TestCase):

    def setUp(self) -> None:
        mod_df = pd.read_csv('../data/MOD17/modeled_results.csv')
        obs_df = pd.read_csv('../data/MOD17/template_example.csv')
        group_df = pd.read_csv('../data/MOD17/group_template.csv')

        # Join all dataframes and filter missing observations
        df = du.merge_model_with_obs(mod_df, obs_df)
        df = du.join_with_groups(df, group_df)
        df = du.filter_nan_obs(df)

        # Nested dict defining parameter space for MCMC calibration process
        params = {'LUE_max': {'min': 1.0, 'max': 3.0, 'guess': 2.0},
                  'VPD_max': {'min': 3000, 'max': 6000, 'guess': 4000},
                  'VPD_min': {'min': 325, 'max': 1500, 'guess': 1000},
                  'T_max': {'min': 0, 'max': 40, 'guess': 20},
                  'T_min': {'min': -10, 'max': 0, 'guess': -5}}

        # Misc arguments used in model calibration
        args = {'nChains': 5,
                'nRuns': 100,
                'nFolds': 2,
                'shuffle': True,
                'parameters': params,
                'filter': 'parVPD_max > parVPD_min & parT_max > parT_min'}

        self.out_dir = './tmp'
        # Create instance of model calibration class
        self.calibrator = CalibrateParameters(df=df,
                                              model=MOD17,
                                              out_dir=self.out_dir,
                                              **args)

        self.calibrator.calibrate()

    def test_calibration_iteration(self):

        self.assertEqual(len(glob.glob(os.path.join(self.out_dir, '*.csv'))), 28)

    def test_calibration_error(self):

        df = pd.read_csv(glob.glob(os.path.join(self.out_dir, '*_training.csv'))[0])
        self.assertNotEqual(np.std(df.like1), 0.0)

    def tearDown(self) -> None:

        for f in glob.glob(os.path.join(self.out_dir, '*.csv')):
            if os.path.isfile(f):
                os.remove(f)


if __name__ == '__main__':
    unittest.main()
