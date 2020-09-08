import pandas as pd


def change_timestamp(df: pd.DataFrame) -> pd.DataFrame:
    """
    :param df: pd.DataFrame with a 'system:time_start' column which contains the date in milliseconds since UNIX epoch.
    :return: pd.DataFrame that contains a 'date' column in YYYY-MM-DD format.
    """
    if 'date' not in df.columns:
        df['date'] = pd.to_datetime(df['system:time_start'], unit='ms').astype('str')
    return df


def merge_model_with_obs(df: pd.DataFrame, obs: pd.DataFrame) -> pd.DataFrame:
    """
    Merge GoogleEarthEngine extracted data with flux tower observations
    :param df: pd.DataFrame containing results from either GeeMod16.py or GeeMod17.py
    :param obs: pd.DataFrame with the corresponding flux tower observations
    :return:
    """
    obs = change_timestamp(obs)
    df = change_timestamp(df)
    df = df.merge(obs, how='left', left_on=['date', 'name'], right_on=['date', 'name'])

    return df


def filter_nan_obs(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filters out missing observations so that model is only calibrated against valid data.
    :param df: pd.DataFrame containing tower observations
    :return: pd.DataFrame with NaN tower data removed.
    """
    df = df[df['target'].notnull()]
    return df


def join_with_groups(df: pd.DataFrame, grp_df: pd.DataFrame = None) -> pd.DataFrame:
    """
    Join a dataframe containing model or observation data with a dataframe that contains PFT grouping info.
    :param df: pd.DataFrame with model or observation data
    :param grp_df: pd.DataFrame with grouping data
    :return: Joined pd.DataFrame between model/observations and groups
    """
    if grp_df is None:
        df['group'] = 'NoGroup'
    else:
        df = df.merge(grp_df, how='left', left_on='name', right_on='name')
    return df
