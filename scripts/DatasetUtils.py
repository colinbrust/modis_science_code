import pandas as pd


def change_timestamp(df: pd.DataFrame) -> pd.DataFrame:
    if 'date' not in df.columns:
        df['date'] = pd.to_datetime(df['system:time_start'], unit='ms').astype('str')
    return df


def merge_model_with_obs(df: pd.DataFrame, obs: pd.DataFrame) -> pd.DataFrame:
    obs = change_timestamp(obs)
    df = change_timestamp(df)
    df = df.merge(obs, how='left', left_on=['date', 'name'], right_on=['date', 'name'])

    return df


def filter_nan_obs(df: pd.DataFrame) -> pd.DataFrame:
    df = df[df['target'].notnull()]
    return df


def join_with_groups(df: pd.DataFrame, grp_df: pd.DataFrame) -> pd.DataFrame:
    df = df.merge(grp_df, how='left', left_on='name', right_on='name')
    return df
