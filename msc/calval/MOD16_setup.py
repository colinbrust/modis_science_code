import spotpy
import numpy as np
import pandas as pd
import argparse
import os
from msc.models.LocalMod16_testing_changes import MOD16
from msc.utils import DatasetUtils as du
from sklearn.model_selection import KFold

BPLUT = {'tmin_open': {'min': 8, 'max': 12, 'guess': 10},
         'tmin_close': {'min': -8.0, 'max': -6.0, 'guess': -7.0},
         'vpd_open': {'min': 0, 'max': 2000, 'guess': 1000},
         'vpd_close': {'min': 1450, 'max': 7000, 'guess': 3000},
         'gl_sh': {'min': 0.0, 'max': 0.12, 'guess': 0.1},
         'gl_e_wv': {'min': 0.0, 'max': 0.12, 'guess': 0.1},
         'Cl': {'min': 0.0013, 'max': 0.012, 'guess': 0.005},
         'rbl_min': {'min': 10, 'max': 110, 'guess': 70},
         'rbl_max': {'min': 20, 'max': 150, 'guess': 120}}


class model(object):

    def __init__(self, tmp):

        self.tmp = tmp
        self.observations = tmp.target.values

    def run(self, tmin_open, tmin_close, vpd_close, vpd_open, gl_sh, gl_e_wv, Cl, rbl_min, rbl_max, beta=None,
            sm_open=None, sm_close=None):

        if tmin_open < tmin_close or vpd_open > vpd_close or rbl_min > rbl_max:
            print('min parameter is greater than max, skipping this parameterization')
            return self.observations * -np.inf

        elif (sm_open and sm_close) and (sm_open < sm_close):
            print('min parameter is greater than max, skipping this parameterization')
            return self.observations * -np.inf

        else:
            out = MOD16(self.tmp, tmin_open=tmin_open, tmin_close=tmin_close, vpd_open=vpd_open, vpd_close=vpd_close,
                        gl_sh=gl_sh, gl_e_wv=gl_e_wv, Cl=Cl, rbl_min=rbl_min, rbl_max=rbl_max, beta=beta,
                        sm_close=sm_close, sm_open=sm_open)

            return out.simulation.values


class spotpy_setup(object):
    def __init__(self, tmp, param_bounds, use_sm):
        self.tmp = tmp
        self.param_bounds = param_bounds
        self.use_sm = use_sm

        if self.use_sm:
            self.param_bounds['sm_open'] = {'min': 0.1, 'max': 1.0, 'guess': 0.5}
            self.param_bounds['sm_close'] = {'min': 0.0, 'max': 0.2, 'guess': 0.05}
        else:
            self.param_bounds['beta'] = {'min': 0, 'max': 1000, 'guess': 250}

        self.model = model(self.tmp)
        self.params = [spotpy.parameter.Uniform(x,
                                                low=self.param_bounds[x]['min'],
                                                high=self.param_bounds[x]['max'],
                                                optguess=self.param_bounds[x]['guess']) for x in self.param_bounds]

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):

        if self.use_sm:
            sm_open = vector[9]
            sm_close = vector[10]
            beta = None

        else:
            sm_open = None
            sm_close = None
            beta = vector[9]

        return self.model.run(tmin_open=vector[0], tmin_close=vector[1], vpd_open=vector[2], vpd_close=vector[3],
                              gl_sh=vector[4], gl_e_wv=vector[5], Cl=vector[6], rbl_min=vector[7], rbl_max=vector[8],
                              beta=beta, sm_open=sm_open, sm_close=sm_close)

    def evaluation(self):
        return self.model.observations

    def objectivefunction(self, simulation, evaluation):
        objectivefunction = -spotpy.objectivefunctions.rmse(evaluation, simulation)
        return objectivefunction


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-mod', '--mod_df', type=str, help='DF containing intermedite model variables')
    parser.add_argument('-obs', '--obs_df', type=str, help='DF containing Ameriflux Observations')
    parser.add_argument('-grp', '--group_df', type=str, help='DF containing PFT groupings for towers')
    parser.add_argument('-out', '--out_dir', type=str, help='Directory to write results to')
    parser.add_argument('-nr', '--n_runs', type=int, help='Number of calibration runs')
    parser.add_argument('--use_sm', dest='use_sm', action='store_true')
    parser.add_argument('--no-use_sm', dest='use_sm', action='store_false')
    parser.add_argument('--mcmc', dest='mcmc', action='store_true')
    parser.add_argument('--no-mcmc', dest='mcmc', action='store_false')
    parser.set_defaults(n_runs=5000)
    parser.set_defaults(use_sm=True)
    parser.set_defaults(mcmc=True)

    parser_args = parser.parse_args()

    # Read in observation, necessary GEE model outputs, and file containing station groupings
    mod_df = pd.read_csv(parser_args.mod_df)
    obs_df = pd.read_csv(parser_args.obs_df)
    group_df = pd.read_csv(parser_args.group_df)

    # Join all data frames and filter missing observations
    df = du.merge_model_with_obs(mod_df, obs_df)
    if 'group' in df.columns:
        df = df.drop(columns=['group'])
    df = du.join_with_groups(df, group_df)
    df = du.filter_nan_obs(df)

    grp_list = np.unique(df['group'])

    for grp in grp_list:
        tmp = df[df['group'] == grp]
        kf = KFold(n_splits=10, shuffle=True)

        fold = 0
        for train, test in kf.split(tmp):
            train_df = tmp.reset_index(drop=True).filter(train, axis=0)
            test_df = tmp.reset_index(drop=True).filter(test, axis=0)
            while True:
                try:
                    if parser_args.mcmc:
                        save_name = os.path.join(parser_args.out_dir,
                                                 'nRuns-' + str(parser_args.n_runs) + '_fold-' + str(fold) +
                                                 '_group-' + grp + '_method-mcmc')
                        sampler = spotpy.algorithms.mcmc(spotpy_setup(tmp=train_df, param_bounds=BPLUT,
                                                                      use_sm=parser_args.use_sm),
                                                         dbname=save_name + '_training', dbformat='csv', save_sim=False)
                        sampler.sample(repetitions=parser_args.n_runs)

                    else:
                        save_name = os.path.join(parser_args.out_dir,
                                                 'nRuns-' + str(parser_args.n_runs) + '_fold-' + str(fold) +
                                                 '_group-' + grp + '_method-demc')
                        sampler = spotpy.algorithms.demcz(spotpy_setup(tmp=train_df, param_bounds=BPLUT,
                                                                       use_sm=parser_args.use_sm),
                                                          dbname=save_name + '_training', dbformat='csv',
                                                          save_sim=False)
                        sampler.sample(repetitions=parser_args.n_runs)

                    params = pd.read_csv(save_name + '_training.csv')
                    params = params.groupby('chain').apply(lambda x: x[x['like1'] == x['like1'].max()])
                    params = params.drop(columns=['chain', 'like1']).reset_index().drop(columns=['level_1']).drop(
                        columns=['chain'])
                    params = params.agg('mean').to_frame().reset_index()
                    params['index'] = params['index'].str.replace('par', '')
                    params = params.set_index('index').to_dict()[0]

                    val = MOD16(test_df, **params)
                    val.to_csv(save_name + '_holdout.csv', index=False)
                    fold += 1

                except KeyError as e:
                    print(e)
                    print('Rerunning with new params')
                    continue
                break
