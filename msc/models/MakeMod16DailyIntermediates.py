from math import exp
import numpy as np
import pandas as pd
import os


def calc_tavg(df):
    df['Tavg'] = (df.tmmn + df.tmmx) / 2 - 273.15
    return df


def calc_saturated_vapor_pressure(temp):
    return 610.7 * exp(17.38 * temp / (temp + 239))


def pa_fun(elevation):
    t1 = 1 - (0.0065 * elevation) / 288.15
    t2 = 9.80665 / (0.0065 * 8.3143 / 0.0289644)
    return 101325.0 * (t1 ** t2)


def set_initial_values(df):
    df = calc_tavg(df)
    tannual = df.groupby('name').apply(lambda x: np.mean(x['Tavg']))
    tannual = pd.DataFrame({'Tannual': tannual}).reset_index()
    df = df.merge(tannual, left_on='name', right_on='name')
    df['albedo'] = df.albedo * 0.001
    df['vpd'] = df.vpd * 1000
    df['ta_day'] = df.Tavg
    df['ta_night'] = df.tmmn - 273.15
    df['rh'] = (df.rmax + df.rmin) / 2
    df['rhmin'] = df.rmin
    df['nightlength'] = 86400 - df.dayl
    df['svp_night'] = df['ta_night'].apply(calc_saturated_vapor_pressure)
    df['svp_day'] = df['ta_day'].apply(calc_saturated_vapor_pressure)
    df['vpd_day'] = df.vpd
    df['vpd_night'] = df.svp_night * (1 - (df.rmin * 0.01))
    df['pa'] = df['elevation'].apply(pa_fun)
    df['swrad'] = df.srad
    df['daylength'] = df.dayl
    df['surf_min'] = df.groupby('name')['sm-surface'].transform('min')
    df['surf_max'] = df.groupby('name')['sm-surface'].transform('max')
    df['rz_min'] = df.groupby('name')['sm-rootzone'].transform('min')
    df['rz_max'] = df.groupby('name')['sm-rootzone'].transform('max')

    return df[['name', 'system:time_start', 'swrad', 'albedo', 'Tavg', 'ta_day', 'ta_night', 'rh', 'rhmin',
               'svp_day', 'svp_night', 'vpd_day', 'vpd_night', 'daylength', 'nightlength', 'pa', 'Fc', 'LAI',
               'Tannual', 'sm-surface', 'sm-rootzone', 'surf_min', 'surf_max', 'rz_min', 'rz_max']]


def calc_rn(temp, albedo, swrad, daylength, daytime=True):
    sbc = 5.67 / 1e8
    t = (273.15 + temp) ** 4
    ea = 1 - 0.261 * exp((-7.77) * 0.0001 * temp * temp)
    rn = ((1 - albedo) * swrad + (ea - 0.97) * t * sbc) * daylength

    if daytime:
        rn = np.where(rn < 0, 0, rn)

    return rn


def calc_soil_heat_flux_day(df):
    df['G_tmp'] = (4.73 * df.ta_day - 20.87) * df.daylength
    df['G_day'] = np.where((df['Tannual'] < 25) & (-8 <= df['Tannual']) & ((df['ta_day'] - df['ta_night']) >= 5),
                           df['G_tmp'], -9999)
    df['G_day'] = np.where(abs(df['G_tmp']) > 0.39 * abs(df['Rn_day']), 0.39 * df['Rn_day'], df['G_day'])
    df['G_day'] = np.where(df['G_day'] == -9999, 0, df['G_day'])
    df['G_day'] = np.where(df['Rn_day'] - df['G_day'] < 0, df['Rn_day'], df['G_day'])
    df = df.drop(columns=['G_tmp'])
    return df


def calc_soil_heat_flux_night(df):
    df['G_tmp'] = (4.73 * df.ta_night - 20.87) * df.nightlength
    df['G_night'] = np.where((df['Tannual'] < 25) & (-8 <= df['Tannual']) & ((df['ta_day'] - df['ta_night']) >= 5),
                             df['G_tmp'], -9999)
    df['G_night'] = np.where(abs(df['G_tmp']) > 0.39 * abs(df['Rn_night']), 0.39 * df['Rn_night'], df['G_night'])
    df['G_night'] = np.where(df['G_night'] == -9999, 0, df['G_night'])
    df['G_night'] = np.where((df['Rn_night'] - df['G_night']) < (-0.5 * df['Rn_day']),
                             df['Rn_night'] + 0.5 * df['Rn_day'], df['G_night'])
    df = df.drop(columns=['G_tmp'])
    return df


def calc_canopy_rad(Rn, Fc):
    return Rn * Fc


def calc_soil_rad(Rn, Fc, G):
    return (1 - Fc) * Rn - G


def calc_wet_soil(rh):
    fwet = (rh * 0.01) ** 4
    fwet = np.where(rh < 70, 0, fwet)
    return fwet


def calc_rho(rh, pa, temp):
    return (0.348444 * pa * 0.01 - rh * (0.00252 * temp - 0.020582)) / (temp + 273.15)


def calc_lhv(temp):
    return (2.501 - 0.002361 * temp) * 1e6


def calc_svp_slope(svp, temp):
    t = temp + 239.0
    return 17.38 * 239.0 * svp / (t * t)


def calc_psychrometric_constant(pa, lamda):
    cp = 1013.0
    epsl = 0.622
    return (cp * pa) / (lamda * epsl)


def calc_rht_resistance(rho, temp):
    t = temp + 273.15
    cp = 1013.0
    sbc = 5.67 / 1e8
    return rho * cp / (4.0 * sbc * t * t * t)


def calc_fsm(sm, smMin, smMax):
    return (sm - smMin) / (smMax - smMin)


def calc_all_variables(out_dir='/mnt/e/Data/GEE_ET/clean/calibration_data/met_drivers',
                       df='/mnt/e/Data/GEE_ET/clean/calibration_data/met_drivers/modeled_results.csv'):

    df = pd.read_csv(df)
    df = set_initial_values(df)
    df['Rn_day'] = df.apply(lambda row: calc_rn(row['ta_day'], row['albedo'], row['swrad'], row['daylength'], True),
                            axis=1)
    df['Rn_night'] = df.apply(lambda row: calc_rn(row['ta_night'], row['albedo'], 0, row['nightlength'], False), axis=1)
    print('Calculated Rnet')
    df = calc_soil_heat_flux_day(df)
    df = calc_soil_heat_flux_night(df)
    print('Calculated soil heat flux')
    df['Ac_day'] = df.apply(lambda row: calc_canopy_rad(row['Rn_day'], row['Fc']), axis=1)
    df['Ac_night'] = df.apply(lambda row: calc_canopy_rad(row['Rn_night'], row['Fc']), axis=1)
    print('Calculated canopy radiation')
    df['Asoil_day'] = df.apply(lambda row: calc_soil_rad(row['Rn_day'], row['Fc'], row['G_day']), axis=1)
    df['Asoil_night'] = df.apply(lambda row: calc_soil_rad(row['Rn_night'], row['Fc'], row['G_night']), axis=1)
    print('Calculated soil radiation')
    df['Fwet_day'] = df.apply(lambda row: calc_wet_soil(row['rh']), axis=1)
    df['Fwet_night'] = df.apply(lambda row: calc_wet_soil(row['rhmin']), axis=1)
    print('Calculated wet surface fraction')
    df['rho_day'] = df.apply(lambda row: calc_rho(row['rh'], row['pa'], row['ta_day']), axis=1)
    df['rho_night'] = df.apply(lambda row: calc_rho(row['rhmin'], row['pa'], row['ta_night']), axis=1)
    print('Calculated rho')
    df['lamda_day'] = df.apply(lambda row: calc_lhv(row['ta_day']), axis=1)
    df['lamda_night'] = df.apply(lambda row: calc_lhv(row['ta_night']), axis=1)
    print('Calculated latent heat of vaporization')
    df['s_day'] = df.apply(lambda row: calc_svp_slope(row['svp_day'], row['ta_day']), axis=1)
    df['s_night'] = df.apply(lambda row: calc_svp_slope(row['svp_night'], row['ta_night']), axis=1)
    print('Calculated svp slope')
    df['gama_day'] = df.apply(lambda row: calc_psychrometric_constant(row['pa'], row['lamda_day']), axis=1)
    df['gama_night'] = df.apply(lambda row: calc_psychrometric_constant(row['pa'], row['lamda_night']), axis=1)
    print('Calculated psychrometric constant')
    df['rrs_day'] = df.apply(lambda row: calc_rht_resistance(row['rho_day'], row['ta_day']), axis=1)
    df['rrs_night'] = df.apply(lambda row: calc_rht_resistance(row['rho_night'], row['ta_night']), axis=1)
    print('Calculated rht resistance')
    df['fSM'] = df.apply(lambda row: calc_fsm(row['sm-surface'], row['surf_min'], row['surf_max']), axis=1)
    df['fSM-rz'] = df.apply(lambda row: calc_fsm(row['sm-rootzone'], row['rz_min'], row['rz_max']), axis=1)
    print('done')
    df['date'] = pd.to_datetime(df['system:time_start'], unit='ms').dt.date.astype('str')
    df = df.drop(columns=['system:time_start', 'sm-surface', 'sm-rootzone'])
    df.to_csv(os.path.join(out_dir, 'modeled_results_intermediates.csv'), index=False)
