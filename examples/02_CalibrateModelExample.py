import pandas as pd
from msc.models.LocalMod17 import MOD17
from msc.models.LocalMod16 import MOD16
from msc.utils import DatasetUtils as du
from msc.calval.CalibrateParameters import CalibrateParameters


# ---------- MOD17 Example ----------
# Read in observation, necessary GEE model outputs, and file containing station groupings
mod_df = pd.read_csv('data/MOD17/modeled_results.csv')
obs_df = pd.read_csv('data/MOD17/template_example.csv')
group_df = pd.read_csv('data/MOD17/group_template.csv')

# Join all dataframes and filter missing observations
df = du.merge_model_with_obs(mod_df, obs_df)
df = du.join_with_groups(df, group_df)
df = du.filter_nan_obs(df)

# Nested dict defining parameter space for MCMC calibration process
params = {'LUE_max': {'min': 1.0, 'max': 3.0, 'guess': 2.0},
          'VPD_max': {'min': 3000, 'max': 6000, 'guess': 4000},
          'VPD_min': {'min': 325, 'max': 1500, 'guess': 1000},
          'T_max':   {'min': 0, 'max': 40, 'guess': 20},
          'T_min':   {'min': -10, 'max': 0, 'guess': -5}}

# Misc arguments used in model calibration
args = {'nChains': 5,
        'nRuns': 5000,
        'nFolds': 5,
        'shuffle': True,
        'parameters': params,
        'filter': 'parVPD_max > parVPD_min & parT_max > parT_min'}

# Create instance of model calibration class
calibrator = CalibrateParameters(df=df,
                                 model=MOD17,
                                 out_dir='data/MOD17/kfold_results',
                                 **args)

# calibrate the model
calibrator.calibrate()

# ---------- MOD16 Example ----------
# # Read in observation, necessary GEE model outputs, and file containing station groupings
# mod_df = pd.read_csv('data/MOD16/modeled_results.csv')
# obs_df = pd.read_csv('data/MOD16/template_example.csv')
# group_df = pd.read_csv('data/MOD16/group_template.csv')
#
# # Join all dataframes and filter missing observations
# df = du.merge_model_with_obs(mod_df, obs_df)
# df = du.join_with_groups(df, group_df)
# df = du.filter_nan_obs(df)
#
# # Nested dict defining parameter space for MCMC calibration process
# params = {'tmin_open':  {'min': 8.0, 'max': 12.5, 'guess': 10},
#           'tmin_close': {'min': -8.0, 'max': -6.0, 'guess': -7.0},
#           'vpd_open':   {'min': 325, 'max': 1500, 'guess': 1000},
#           'vpd_close':  {'min': 1450, 'max': 6500, 'guess': 3000},
#           'gl_sh':      {'min': 0.005, 'max': 0.12, 'guess': 0.1},
#           'gl_e_wv':    {'min': 0.005, 'max': 0.12, 'guess': 0.1},
#           'Cl':         {'min': 0.0013, 'max': 0.013, 'guess': 0.005},
#           'rbl_min':    {'min': 10, 'max': 105, 'guess': 70},
#           'rbl_max':    {'min': 20, 'max': 150, 'guess': 120}}
#
# # Misc arguments used in model calibration
# args = {'nChains': 5,
#         'nRuns': 100,
#         'nFolds': 5,
#         'shuffle': True,
#         'parameters': params,
#         'filter': 'parvpd_close > parvpd_open & parrbl_min < parrbl_max'}
#
# # Create instance of model calibration class
# calibrator = CalibrateParameters(df=df,
#                                  model=pd_MOD16,
#                                  out_dir='data/MOD16/kfold_results',
#                                  **args)
#
# # calibrate the model
# calibrator.calibrate()
