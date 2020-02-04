import pandas as pd
from typing import Optional


class Dataset(object):

    def __init__(self, df: pd.DataFrame, obs_df: Optional[pd.DataFrame] = None) -> None:
        self.df = df
        self.df = self._change_timestamp(self.df)
        if obs_df is not None:
            obs_df = self._change_timestamp(obs_df)
            self.merge_model_with_obs(obs_df)
        self.target = None
        self.simulation = None
        self.fill_target()
        self.fill_simulation()

    @staticmethod
    def _change_timestamp(df):
        if 'date' not in df.columns:
            df['date'] = pd.to_datetime(df['system:time_start'], unit='ms').astype('str')
        return df

    def merge_model_with_obs(self, obs) -> pd.DataFrame:
        obs = self._change_timestamp(obs)
        self.df = self.df.merge(obs, how='left', left_on=['date', 'name'], right_on=['date', 'name'])
        self.fill_simulation()
        self.fill_target()

    def filter_nan_obs(self):
        self.df = self.df[self.df['target'].notnull()]
        self.fill_simulation()
        self.fill_target()

    def fill_simulation(self):
        self.simulation = self.df['simulation'].values if 'simulation' in self.df.columns else None

    def fill_target(self):
        self.target = self.df['target'].values if 'target' in self.df.columns else None

    def join_with_groups(self, grp_df):
        self.df = self.df.merge(grp_df, how='left', left_on='name', right_on='name')

    def __repr__(self):
        names = self.df.name.values.tolist()
        years = self.df.year.values.tolist()

        out = [name + ': ' + str(year) for name, year in zip(names, years)]
        out = out[0:3] + ['  ......  '] + out[-3:]
        out = '\n'.join(out)
        return 'Dataset with data for following sites/years:\n'+out
