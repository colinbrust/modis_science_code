import spotpy
import pandas as pd
import numpy as np
import os
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

    def train_and_test(self, df, group):

        fold = 0
        kf = KFold(n_splits=self.nFolds, shuffle=self.shuffle)
        for train, test in kf.split(df):

            train_df = df.filter(train, axis=0)
            test_df = df.filter(test, axis=0)
            save_name = os.path.join()

            model = Calibration(param_dict=self.param_dict, model=self.model, df=train_df)

            sampler = spotpy.algorithms.demcz(model,
                                              dbname=out_dir + '/' + str(nrun) + '_' + pft + '_fold' + str(fold) + '_' +
                                                     str(nChains), dbformat='csv', save_sim=False, parallel='mpc')
            sampler.sample(repetitions=nrun, nChains=nChains)


            validation = et_calc(test_df, tmin_open=bplut['tmin_open'][pft_dict[pft]],
                                 tmin_close=bplut['tmin_close'][pft_dict[pft]], vpd_open=params['parvpd_open'],
                                 vpd_close=params['parvpd_close'], sm_close=params['parsm_close'],
                                 sm_open=params['parsm_open'], gl_sh=params['pargl_sh'], gl_e_wv=params['pargl_e_wv'],
                                 Cl=params['parCl'], rbl_min=params['parrbl_min'], rbl_max=params['parrbl_max'])

            validation.to_csv(os.path.join(out_dir, pft + '_' + str(fold) + '_holdout.csv'), index=False)

            params['pft'] = pft
            err_df.append(params)
            i += 1

        return pd.DataFrame(err_df)