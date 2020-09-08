import pandas as pd
import numpy as np

"""
To calculate GPP, a pandas dataframe is required that has daily minimum temperature, vpd, solar radiation, and fpar in
columns labelled 'tmmn', 'vpd', 'srad', and 'fpar', respectively. An example dataframe can be found in 
'data/MOD17/modeled_results.csv'
"""


def calc_lue(df, LUE_max, VPD_max, VPD_min, T_max, T_min):
    """
    :param df: pd.DataFrame that contains 'tmmn' and 'vpd' columns that correspond to the daily minimum temperature in K
        and the vapor pressure deficit in kPa, respectively
    :param LUE_max: Maximum light use efficiency parameter
    :param VPD_max: VPD at which GPP is maximized parameter
    :param VPD_min: VPD at which GPP stops parameter
    :param T_max: Temperature at which GPP is maximized parameter
    :param T_min: Temperature at which GPP stops parameter
    :return: pd.DataFrame with linear GPP scalars and scaled light use efficiency.
    """
    df['Tscalar'] = ((df.tmmn - 273.15) - T_min) / (T_max - T_min)
    df['Tscalar'] = np.clip(df.Tscalar, 0, 1)

    df['VPDscalar'] = (VPD_max - (df.vpd * 1000)) / (VPD_max - VPD_min)
    df['VPDscalar'] = np.clip(df.VPDscalar, 0, 1)

    df['lue'] = LUE_max * df.Tscalar * df.VPDscalar

    return df


# Function to calculate GPP
def MOD17(df, LUE_max, VPD_max, VPD_min, T_max, T_min):
    """

    :param df: pd.DataFrame that contains 'tmmn', 'fpar', 'srad', and 'vpd' columns that correspond to the daily minimum
        temperature in K, fraction of vegetation cover, solar radiation in W/m^2, and the vapor pressure deficit in kPa.
    :param LUE_max: Maximum light use efficiency parameter
    :param VPD_max: VPD at which GPP is maximized parameter
    :param VPD_min: VPD at which GPP stops parameter
    :param T_max: Temperature at which GPP is maximized parameter
    :param T_min: Temperature at which GPP stops parameter
    :return: pd.DataFrame with GPP, labelled as 'simulation'
    """
    df = calc_lue(df, LUE_max, VPD_max, VPD_min, T_max, T_min)
    df['par'] = df.srad * 0.45 * 60 * 60 * 24 * 0.000001
    df['simulation'] = df.fpar * df.par * df.lue

    return df
