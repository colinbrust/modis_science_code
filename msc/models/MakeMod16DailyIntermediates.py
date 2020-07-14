from math import exp
import numpy as np
import pandas as pd

def calc_tavg(df):

    df['Tavg'] = (df.tmmn + df.tmmx)/2-273.15
    return df


def calc_saturated_vapor_pressure(temp):
    return 610.7*exp(17.38*temp/(temp+239))


def pa_fun(elevation):
    t1 = 1-(0.0065*elevation)/288.15
    t2 = 9.80665/(0.0065*8.3143/0.0289644)
    return 101325.0 * (t1 ** t2)


def set_initial_values(df):

    df = calc_tavg(df)
    tannual = df.groupby('name').apply(lambda x: np.mean(x['Tavg']))
    tannual = pd.DataFrame({'Tannual':tannual}).reset_index()
    df = df.merge(tannual, left_on='name', right_on='name')
    df['albedo'] = df.albedo * 0.001
    df['vpd'] = df.vpd * 1000
    df['ta_day'] = df.Tavg
    df['ta_night'] = df.tmmn - 273.15
    df['rh_day'] = (df.rmax + df.rmin)/2
    df['rh_night'] = df.rmin
    df['nightlength'] = 86400 - df.dayl
    df['svp_night'] = df['ta_night'].apply(calc_saturated_vapor_pressure)
    df['svp_day'] = df['ta_day'].apply(calc_saturated_vapor_pressure)
    df['vpd_day'] = df.vpd
    df['vpd_night'] = df.svp_night * (1 - (df.rmin * 0.01))
    df['pa'] = df['elevation'].apply(pa_fun)
    df['swrad'] = df.srad
    df['daylength'] = df.dayl

    return df[['name', 'system:time_start', 'swrad', 'albedo', 'Tavg', 'ta_day', 'ta_night', 'rh_day', 'rh_night',
               'svp_day', 'svp_night', 'vpd_day', 'vpd_night', 'daylength', 'nightlength', 'pa', 'Fc', 'LAI']]


def calc_rn(temp, albedo, swrad, daylength, daytime=True):

    sbc = 5.67 / 1e8
    t = (273.15 + temp) ** 4
    ea = 1 - 0.261 * exp((-7.77) * 0.0001 * temp * temp)
    rn = ((1 - albedo) * swrad + (ea - 0.97) * t * sbc) * daylength

    if daytime:
        rn = np.where(rn < 0, 0, rn)

    return rn


def calc_soil_heat_flux(tmin_close, rn, temp, daylength, ta_day, ta_night):

    gsoil = (4.73 * temp - 20.87) * daylength
    if tmin_close <=

def calc_all_variables(df):

    df = set_initial_values(df)
    df['Rn_day'] = df.apply(lambda row: calc_rn(row['ta_day'], row['albedo'], row['swrad'], row['daylength'], True),
                            axis=1)
    df['Rn_night'] = df.apply(lambda row: calc_rn(row['ta_night'], row['albedo'], 0, row['nightlength'], False), axis=1)