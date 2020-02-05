import pandas as pd
from scripts.LocalMod16 import MOD16
import scripts.DatasetUtils as du
import scripts.EvaluateCalibration as ec

# read read in gee extracted data necessary to run model locally
df = pd.read_csv('data/MOD16/modeled_results.csv')

# read dataframe that maps site names to pft groups
group_df = pd.read_csv('data/MOD16/group_template.csv')
df = du.join_with_groups(df, group_df)

# assign directory that has results of calibration
kfold_dir = 'data/MOD16/kfold_results'

# Calculate holdout error from calibration
ec.calc_holdout_error(kfold_dir=kfold_dir, group_col='group')

# Get optimal parameters from calibration process
ec.get_optimal_params(kfold_dir, 'parvpd_close > parvpd_open & parrbl_min < parrbl_max')

# Get dataframe with model results using new parameter values
new_data = ec.run_all_with_opt_params(df=df, kfold_dir=kfold_dir, model=MOD16,
                                      cond_filter='parvpd_close > parvpd_open & parrbl_min < parrbl_max')