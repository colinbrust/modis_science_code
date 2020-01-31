import pandas as pd


class Dataset(object):

    def __init__(self, obs_df: pd.DataFrame, mod_df: pd.DataFrame, target: str, simulation: str) -> None:
        self.obs_df = obs_df
        self.mod_df = mod_df
        self._change_modeled_timestamp()
        self.target = target
        self.simulation = simulation
        self.df = self._merge_model_with_obs()

    def _change_modeled_timestamp(self):
        if 'date' not in self.mod_df.columns:
            self.mod_df['date'] = pd.to_datetime(self.mod_df['system:time_start'], unit='ms').astype('str')

    def _merge_model_with_obs(self) -> pd.DataFrame:
        df = self.mod_df.merge(self.obs_df, how='left', left_on=['date', 'name'], right_on=['date', 'name'])