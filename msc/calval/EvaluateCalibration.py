import pandas as pd
from glob import glob
import os
from spotpy.objectivefunctions import rmse, bias, rsquared
from typing import Callable


def calc_holdout_error(kfold_dir: str, group_col: str = None) -> pd.DataFrame:
    """
    Function to calculate error using holdout data that was not used for model calibration.
    :param kfold_dir: Directory that k-fold cross validation results are put into
    :param group_col: data frame column to group error metrics by. 'name' if you want error by station, 'group' if you
        want error by 'group' column, None if you want error across all sites.
    :return: pd.DataFrame that gives RMSE, bias and R2 for each group
    """
    df = pd.DataFrame()
    for f in glob(os.path.join(kfold_dir, '*_holdout.csv')):
        df = pd.concat([df, pd.read_csv(f)], ignore_index=True)

    if group_col is not None:
        rmse_dict = {'rmse': df.groupby(group_col).apply(lambda x: rmse(x['target'], x['simulation'])).to_dict()}
        bias_dict = {'bias': df.groupby(group_col).apply(lambda x: bias(x['target'], x['simulation'])).to_dict()}
        r2_dict = {'r2': df.groupby(group_col).apply(lambda x: rsquared(x['target'], x['simulation'])).to_dict()}

        return pd.DataFrame({**rmse_dict, **bias_dict, **r2_dict})

    else:
        rmse_dict = {'rmse': rmse(df['target'], df['simulation'])}
        bias_dict = {'bias': bias(df['target'], df['simulation'])}
        r2_dict = {'r2': rsquared(df['target'], df['simulation'])}

        return pd.DataFrame({**rmse_dict, **bias_dict, **r2_dict}, index=[0])


def get_optimal_params(kfold_dir: str, cond_filter: str = None) -> pd.DataFrame:
    """
    Function to find the optimal parameters across all training folds of the data.
    :param kfold_dir:  Directory that k-fold cross validation results are put into
    :param cond_filter: Conditional filter for model parameters.  See examples/02_CalibrateModelExample.py for more
        details.
    :return: pd.DataFrame that contains mean of optimal parameters across all folds for all groups
    """
    df = pd.DataFrame()
    for f in glob(os.path.join(kfold_dir, '*_training.csv')):
        tmp = pd.read_csv(f)
        tmp['fold'] = [x for x in f.split('_') if 'fold-' in x][0].replace('fold-', '')
        tmp['group'] = [x for x in f.split('_') if 'group-' in x][0].replace('group-', '')
        df = pd.concat([df, tmp], ignore_index=True)

    df = df.query(cond_filter) if cond_filter is not None else df
    df = df.groupby(['chain', 'fold', 'group']).apply(lambda x: x[x['like1'] == x['like1'].max()])
    df = df.reset_index(drop=True).groupby('group').agg('mean').reset_index()
    df.columns = [x.replace('par', '') for x in df.columns]

    return df.drop(columns=['like1', 'chain'])


def run_all_with_opt_params(df: pd.DataFrame, kfold_dir: str, model: Callable, cond_filter: str = None, premade: pd.DataFrame = None) -> pd.DataFrame:
    """
    Runs a model using optimal calibrated parameters.
    :param df: pd.DataFrame containing all necessary data to run model.
    :param kfold_dir: Directory that k-fold cross validation results are put into
    :param model: Model to run.
    :param cond_filter: Conditional filter for model parameters.  See examples/02_CalibrateModelExample.py for more
        details.
    :return:
    """
    params = premade if premade is not None else get_optimal_params(kfold_dir=kfold_dir, cond_filter=cond_filter)
    groups = params['group'].values

    if len(groups) == 1:
        opt_params = params[params['group'] == groups[0]].drop(columns='group').to_dict('records')[0]
        return model(df, **opt_params)
    else:
        df_out = pd.DataFrame()
        for group in groups:
            opt_params = params[params['group'] == group].drop(columns='group').to_dict('records')[0]
            tmp = df[df['group'] == group]
            tmp = model(tmp, **opt_params)
            df_out = pd.concat([tmp, df_out], ignore_index=True)
        return df_out
