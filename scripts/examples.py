import ee
import pandas as pd
from scripts.ee_MOD16 import MOD16
from scripts.ExtractGeeData import ExtractGeeData
from scripts.local_MOD16 import MOD16 as local_MOD16
from scripts.Dataset import Dataset
from scripts.CalibrateParameters import CalibrateParameters

# --------------------------Calibrate MOD16 Example--------------------------------

mod_df = pd.read_csv('data/modeled_results.csv')
obs_df = pd.read_csv('data/template_example.csv')
group_df = pd.read_csv('data/group_template.csv')

dat = Dataset(df=mod_df, obs_df=obs_df)
dat.join_with_groups(group_df)
dat.filter_nan_obs()

params = {'tmin_open':  {'min': 8.0, 'max': 12.5, 'guess': 10},
          'tmin_close': {'min': -8.0, 'max': -6.0, 'guess': -7.0},
          'vpd_open':   {'min': 325, 'max': 1500, 'guess': 1000},
          'vpd_close':  {'min': 1450, 'max': 6500, 'guess': 3000},
          'gl_sh':      {'min': 0.005, 'max': 0.12, 'guess': 0.1},
          'gl_e_wv':    {'min': 0.005, 'max': 0.12, 'guess': 0.1},
          'Cl':         {'min': 0.0013, 'max': 0.013, 'guess': 0.005},
          'rbl_min':    {'min': 10, 'max': 105, 'guess': 70},
          'rbl_max':    {'min': 20, 'max': 150, 'guess': 120}}

args = {'nChains': 5,
        'nRuns': 100,
        'nFolds': 10,
        'shuffle': True,
        'parameters': params}

calibrator = CalibrateParameters(dataset=dat,
                                 model=local_MOD16,
                                 out_dir='/mnt/e/PycharmProjects/modis_science_code/data',
                                 **args)

calibrator.calibrate()

# # ----------------------------MOD16 in EarthEngine Example-------------------------------------
#
# year = 2015
# meteo = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')
#
# # Here the Daymet daylength product is used to differentiate between daytime and nighttime ET.
# # However, Daymet is only available for North America. In the original MOD16 algorithm, daytime
# # is defined as the period of time when incoming SW radiation is > 10 W/m^2.
# dayl = ee.ImageCollection('NASA/ORNL/DAYMET_V3').select('dayl')
# elev = ee.Image('USGS/NED')
#
# extractor = ExtractGeeData(template='../data/template_example.csv',
#                            model=MOD16,
#                            out_dir='../data',
#                            meteorology=meteo,
#                            daylength=dayl,
#                            elev=elev)
#
# extractor.run_model(restart=False)
