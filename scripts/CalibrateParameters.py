import spotpy
import pandas as pd

class CalibrateP

    def split_data(df, folds):

        split_len = int(len(df) / folds)
        remaining = len(df) - split_len * folds

        base = [x for x in range(folds)] * split_len
        base.sort()

        splits = base + [base[-1]] * remaining
        df['splits'] = splits

        for i in range(folds):
            train = df[df.splits != i]
            test = df[df.splits == i]

            yield train, test

        # if len(sites) == 1:
        #     split_len = round(len(df) / folds) - 1
        #     remaining = len(df) - split_len * folds
        #
        #     splits = ([0] * split_len) + ([1] * split_len) + ([2] * split_len) + ([2] * remaining)
        #
        #     df['splits'] = splits
        #
        #     for i in range(folds):
        #         train = df[df.splits != i]
        #         test = df[df.splits == i]
        #
        #         yield train, test
        #
        # else:
        #
        #     for site in sites:
        #         train = df[df.name != site]
        #         test = df[df.name == site]
        #
        #         yield train, test

    def train_and_test(out_dir='/home/colin.brust/workspace/calibration/kfold_results', nrun=100, nChains=10,
                       f_dir='/home/colin.brust/workspace/data/GEE_ET/pft_calibration', pft='C3'):

        pft_dict = {'CRO': 0, 'ENF': 1, 'EBF': 2, 'DNF': 3, 'DBF': 4,
                    'MF': 5, 'CSH': 6, 'OSH': 7, 'WSA': 8, 'SAV': 9, 'GRA': 10, 'C3': 0, 'C4': 0}

        bplut = {'tmin_open': [12.02, 8.31, 9.09, 10.44, 9.94, 9.50, 8.61, 8.80, 11.39, 11.39, 12.02],
                 'tmin_close': [-8.0, -8.0, -8.0, -8.0, -6.0, -7.0, -8.0, -8.0, -8.0, -8.0, -8.0]}

        df = pd.read_csv(f_dir + '/' + pft + '.csv')

        err_df = []

        i = 0
        for train, test in split_data(df=df, folds=10):
            fold = i
            # train_df = df.filter(train, axis=0)
            # test_df = df.filter(test, axis=0)

            train_df = train
            test_df = test

            model = Calibration(f_dir=f_dir, pft=pft, df=train_df)
            sampler = spotpy.algorithms.demcz(model,
                                              dbname=out_dir + '/' + str(nrun) + '_' + pft + '_fold' + str(fold) + '_' +
                                                     str(nChains), dbformat='csv', save_sim=False, parallel='mpc')
            sampler.sample(repetitions=nrun, nChains=nChains)

            params = run_best_params(
                param_file=out_dir + '/' + str(nrun) + '_' + pft + '_fold' + str(fold) + '_' + str(nChains) + '.csv')

            validation = et_calc(test_df, tmin_open=bplut['tmin_open'][pft_dict[pft]],
                                 tmin_close=bplut['tmin_close'][pft_dict[pft]], vpd_open=params['parvpd_open'],
                                 vpd_close=params['parvpd_close'], sm_close=params['parsm_close'],
                                 sm_open=params['parsm_open'], gl_sh=params['pargl_sh'], gl_e_wv=params['pargl_e_wv'],
                                 Cl=params['parCl'], rbl_min=params['parrbl_min'], rbl_max=params['parrbl_max'])

            validation.to_csv(os.path.join(out_dir, pft + '_' + str(fold) + '_holdout.csv'), index=False)

            params['rmse'] = spotpy.objectivefunctions.rmse(validation.ET_flux.values, validation.ET_mod.values)
            params['corr'] = spotpy.objectivefunctions.rsquared(validation.ET_flux.values, validation.ET_mod.values)
            params['bias'] = spotpy.objectivefunctions.bias(validation.ET_flux.values, validation.ET_mod.values)
            params['rrmse'] = spotpy.objectivefunctions.rrmse(validation.ET_flux.values, validation.ET_mod.values)
            params['pft'] = pft
            err_df.append(params)
            i += 1

        return pd.DataFrame(err_df)