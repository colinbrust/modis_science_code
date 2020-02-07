import pandas as pd
import numpy as np


def calc_lue(df, LUE_max, VPD_max, VPD_min, T_max, T_min):

    df['Tscalar'] = ((df.tmmn - 273.15) - T_min) / (T_max - T_min)
    df['Tscalar'] = np.clip(df.Tscalar, 0, 1)

    df['VPDscalar'] = (VPD_max - (df.vpd * 1000)) / (VPD_max - VPD_min)
    df['VPDscalar'] = np.clip(df.VPDscalar, 0, 1)

    df['lue'] = LUE_max * df.Tscalar * df.VPDscalar

    return df


# Function to calculate GPP
def MOD17(df, LUE_max, VPD_max, VPD_min, T_max, T_min):

    df = calc_lue(df, LUE_max, VPD_max, VPD_min, T_max, T_min)
    df['par'] = df.srad * 0.45 * 60 * 60 * 24 * 0.000001
    df['simulation'] = df.fpar * df.par * df.lue

    return df
