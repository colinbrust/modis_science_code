import spotpy
import numpy as np
import os
import pandas as pd
from sklearn.model_selection import KFold
from typing import Callable
from scripts.Dataset import Dataset
from scripts.Calibration import Calibration


class CalibrateParameters(object):

    def __init__(self, dataset: Dataset, model: Callable, out_dir: str, **kwargs) -> None:
        self.dataset = dataset
        self.model = model
        self.out_dir = out_dir
        self.kwargs = kwargs
        self.nChains = kwargs['nChains'] if 'nChains' in kwargs else 10
        self.nRuns = kwargs['nRuns'] if 'nRuns' in kwargs else 5000
        self.nFolds = kwargs['nFolds'] if 'nFolds' in kwargs else 10
        self.shuffle = kwargs['shuffle'] if 'shuffle' in kwargs else True
        self.grp_list = np.unique(self.dataset.df['group']) if 'group' in self.dataset.df.columns else None

        if 'parameters' not in kwargs:
            raise KeyError("parameters not in kwargs. Parameters must be present for calibration")
        self.param_dict = kwargs['parameters']

    def _find_opt_params(self, param_file):

        params = pd.read_csv(param_file)
        if 'filter' in self.kwargs['filter']:
            params = params.query(self.kwargs['filter'])

        params = params.groupby('chain').apply(lambda x: x[x['like1'] == x['like1'].max()])
        params = params.drop(columns=['chain', 'like1']).reset_index().drop(columns=['level_1']).drop(columns=['chain'])
        params = params.agg('mean').to_frame().reset_index()
        params['index'] = params['index'].str.replace('par', '')

        return params.set_index('index').to_dict()[0]

    def _train_and_test(self, df, group):

        if group is not None:
            df = df[df['group'] == group]

        fold = 0
        kf = KFold(n_splits=self.nFolds, shuffle=self.shuffle)
        for train, test in kf.split(df):

            train_df = Dataset(df.filter(train, axis=0))
            test_df = Dataset(df.filter(test, axis=0))
            grp_name = group if group is not None else 'NoGroup'
            save_name = os.path.join(self.out_dir, 'nRuns-'+str(self.nRuns)+'_fold-'+str(fold)+'_group-'+grp_name)

            model = Calibration(param_dict=self.param_dict, model=self.model, dataset=train_df)
            sampler = spotpy.algorithms.demcz(model, dbname=save_name, dbformat='csv', save_sim=False)
            sampler.sample(repetitions=self.nRuns, nChains=self.nChains)

            params = self._find_opt_params(save_name+'.csv')
            validation = self.model(test_df, **params)
            validation = validation.df
            validation.to_csv(os.path.join(self.out_dir, 'fold-'+str(fold)+'_group-'+grp_name+'_holdout.csv'),
                              index=False)

            fold += 1

    def calibrate(self):

        if self.grp_list is None:
            self._train_and_test(self.dataset.df.copy(), None)
        else:
            for grp in self.grp_list:
                self._train_and_test(self.dataset.df.copy(), grp)
