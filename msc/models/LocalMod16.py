import numpy as np

SBC = 5.67 / 1e8  # (W/(m2 K4)) Stefan-Boltzmann constant
Cp = 1013.0
epsl = 0.622


# calculate latent heat flux for wet canopy
def LEwetfun_day(df, gl_sh, gl_e_wv):

    df['rhc'] = (1.0/(gl_sh*df.LAI*df.Fwet_day))
    df['t'] = df.ta_day + 273.15
    df['rrc'] = df.rho_day*Cp/(4.0*SBC*df.t*df.t*df.t)
    df['rhrc'] = df.rhc*df.rrc/(df.rhc+df.rrc)
    df['rvc'] = 1.0/(gl_e_wv*df.LAI*df.Fwet_day)
    df['LEwet_day'] = (df.s_day * df.Ac_day + df.rho_day * Cp * df.vpd_day * df.Fc * df.daylength / df.rhrc) * df.Fwet_day /\
                      (df.s_day + (df.pa * Cp * df.rvc / (df.lamda_day * epsl * df.rhrc)))

    df = df.drop(columns=['t', 'rhc', 'rrc', 'rhrc', 'rvc'])
    return df


def LEwetfun_night(df, gl_sh, gl_e_wv):

    df['rhc'] = 1.0/(gl_sh*df.LAI*df.Fwet_night)
    df['t'] = df.ta_night + 273.15
    df['rrc'] = df.rho_night*Cp/(4.0*SBC*df.t*df.t*df.t)
    df['rhrc'] = df.rhc*df.rrc/(df.rhc+df.rrc)
    df['rvc'] = 1.0/(gl_e_wv*df.LAI*df.Fwet_night)
    df['LEwet_night'] = (df.s_night * df.Ac_night + df.rho_night * Cp * df.vpd_night * df.Fc * df.nightlength / df.rhrc)\
                         * df.Fwet_night / (df.s_night + (df.pa * Cp * df.rvc / (df.lamda_night * epsl * df.rhrc)))

    df = df.drop(columns=['t', 'rhc', 'rrc', 'rhrc', 'rvc'])
    return df


# ------------------------------------------------Plant transpiration-----------------------
# m(vpd)
def mVPD_day(df, vpd_close, vpd_open):

    df['mvpd'] = (vpd_close-df.vpd_day)/(vpd_close-vpd_open)
    df['mvpd'] = np.clip(df.mvpd, 0, 1)

    return df


def mTemp_day(df, tmin_close, tmin_open):

    df['mtemp'] = (df.ta_day - tmin_close)/(tmin_open-tmin_close)
    df['mtemp'] = np.clip(df.mtemp, 0, 1)

    return df


def mSM(df, sm_open, sm_close):

    df['msm'] = (df['sm-rootzone'] - sm_close)/(sm_open - sm_close)
    df['msm'] = np.clip(df.msm, 0, 1)

    return df


# daytime plant transpiration
def LEtrans_day(df, Cl, gl_sh, vpd_open, vpd_close, tmin_open, tmin_close, sm_open, sm_close):

    df['ta'] = ((df.ta_day+273.15)/293.15)**1.75
    dat = df
    dat['rcorr'] = 1.0/((101300/dat.pa)*dat.ta)

    dat = mVPD_day(dat, vpd_close=vpd_close, vpd_open=vpd_open)
    dat = mTemp_day(dat, tmin_close=tmin_close, tmin_open=tmin_open)

    if sm_open is not None and sm_close is not None:
        dat = mSM(df, sm_open=sm_open, sm_close=sm_close)
        dat['Gs1'] = Cl * dat.mvpd * dat.mtemp * dat.msm * dat.rcorr

    else:
        dat['Gs1'] = Cl*dat.mvpd*dat.mtemp*dat.rcorr
    dat['Gcu'] = dat.rcorr*0.00001
    dat['Gs2'] = gl_sh
    dat['Cc'] = dat.LAI*(1.0-dat.Fwet_day)*(dat.Gs2*(dat.Gs1+dat.Gcu))/(dat.Gs1+dat.Gs2+dat.Gcu)

    dat['Cc'] = np.where((dat['LAI'] > 0) & ((1.0 - dat['Fwet_day']) > 0), dat['Cc'], 0)
    dat['rs'] = 1.0/dat.Cc
    dat['rh'] = 1.0/gl_sh
    dat['rr'] = dat.rho_day*Cp/(4.0*SBC*(dat.ta_day+273.15)*(dat.ta_day+273.15)*(dat.ta_day+273.15))
    dat['ra'] = (dat.rh*dat.rr)/(dat.rh+dat.rr)
    dat['LEtrans'] = (dat.s_day*dat.Ac_day+dat.rho_day*Cp*dat.vpd_day*dat.Fc*dat.daylength/dat.ra)*(1.-dat.Fwet_day)/(dat.s_day+dat.gama_day*(1+dat.rs/dat.ra))

    df['LEtrans_day'] = dat.LEtrans
    return df


# nighttime plant transpiration
def LEtrans_night(df, gl_sh):

    df['ta'] = ((df.ta_night+273.15)/293.15)**1.75
    dat = df

    dat['rcorr'] = 1.0/((101300/dat.pa)*dat.ta)

    dat['Gs1'] = 0
    dat['Gcu'] = dat.rcorr*0.00001
    dat['Gs2'] = gl_sh
    dat['Cc'] = dat.LAI*(1.0-dat.Fwet_night)*(dat.Gs2*(dat.Gs1+dat.Gcu))/(dat.Gs1+dat.Gs2+dat.Gcu)

    dat['Cc'] = np.where((dat['LAI'] > 0) & ((1.0 - dat['Fwet_night']) > 0), dat['Cc'], 0)
    dat['rs'] = 1.0/dat.Cc
    dat['rh'] = 1.0/gl_sh
    dat['rr'] = dat.rho_night*Cp/(4.0*SBC*(dat.ta_night+273.15)*(dat.ta_night+273.15)*(dat.ta_night+273.15))
    dat['ra'] = (dat.rh*dat.rr)/(dat.rh+dat.rr)
    dat['LEtrans'] = (dat.s_night*dat.Ac_night+dat.rho_night*Cp*dat.vpd_night*dat.Fc*dat.nightlength/dat.ra)*(1.-dat.Fwet_night)/(dat.s_night+dat.gama_night*(1+dat.rs/dat.ra))

    df['LEtrans_night'] = dat.LEtrans
    return df


# potential plant transpiration is calculated following the Priestley-Taylor
def PLE_plant_day(df):
    df['PLE_plant_day'] = (1.26*df.s_day*df.Ac_day*(1-df.Fwet_day)/(df.s_day+df.gama_day))
    return df

def PLE_plant_night(df):
    df['PLE_plant_night'] = (1.26 * df.s_night * df.Ac_night * (1 - df.Fwet_night) / (df.s_night + df.gama_night))
    return df


# ------------------------------------------Soil surface evaporation calculations--------------
def rtot_day(df, rbl_max, rbl_min, vpd_close, vpd_open):

    rblmax = rbl_max
    rblmin = rbl_min
    vpdclose = vpd_close
    vpdopen = vpd_open

    df['ta'] = ((df.ta_day + 273.15) / 293.15) ** 1.75
    dat = df
    dat['rcorr'] = 1.0 / ((101300 / dat.pa) * dat.ta)
    dat['rtotc'] = rblmax-((rblmax-rblmin)*(vpdclose-dat.vpd_day))/(vpdclose-vpdopen)

    dat['rtotc'] = np.where((dat['vpd_day'] <= vpdopen), rblmax, dat['rtotc'])
    dat['rtotc'] = np.where((dat['vpd_day'] >= vpdclose), rblmin, dat['rtotc'])
    dat['rtot_day'] = dat.rtotc*dat.rcorr
    df['rtot_day'] = dat.rtot_day

    return df


def rtot_night(df, rbl_max, rbl_min, vpd_close, vpd_open):

    rblmax = rbl_max
    rblmin = rbl_min
    vpdclose = vpd_close
    vpdopen = vpd_open

    df['ta'] = ((df.ta_night + 273.15) / 293.15) ** 1.75
    dat = df
    dat['rcorr'] = 1.0 / ((101300 / dat.pa) * dat.ta)
    dat['rtotc'] = rblmax-((rblmax-rblmin)*(vpdclose-dat.vpd_night))/(vpdclose-vpdopen)

    dat['rtotc'] = np.where((dat['vpd_night'] <= vpdopen), rblmax, dat['rtotc'])
    dat['rtotc'] = np.where((dat['vpd_night'] >= vpdclose), rblmin, dat['rtotc'])
    dat['rtot_night'] = dat.rtotc*dat.rcorr
    df['rtot_night'] = dat.rtot_night

    return df


def LEwetsoil_day(df):

    df['ras'] = df.rtot_day*df.rrs_day/(df.rtot_day+df.rrs_day)
    df['LEwetsoil_day'] = (df.s_day*df.Asoil_day+df.rho_day*Cp*(1.0-df.Fc)*df.vpd_day*df.daylength/df.ras)*df.Fwet_day / (df.s_day+df.gama_day*df.rtot_day/df.ras)

    return df

def LEwetsoil_night(df):

    df['ras'] = df.rtot_night*df.rrs_night/(df.rtot_night+df.rrs_night)
    df['LEwetsoil_night'] = (df.s_night*df.Asoil_night+df.rho_night*Cp*(1.0-df.Fc)*df.vpd_night*df.nightlength/df.ras)*df.Fwet_night / (df.s_night+df.gama_night*df.rtot_night/df.ras)

    return df


# potential soil evaporation
def PLE_soil_night(df):

    df['ras']=df.rtot_night * df.rrs_night / (df.rtot_night + df.rrs_night)
    df['PLEsoil_night'] = (df.s_night*df.Asoil_night+df.rho_night*Cp*(1.0-df.Fc)*df.vpd_night*df.nightlength/df.ras) * \
                          (1.0-df.Fwet_night)/(df.s_night+df.gama_night*df.rtot_night/df.ras)

    return df


def PLE_soil_day(df):

    df['ras'] = df.rtot_day * df.rrs_day / (df.rtot_day + df.rrs_day)
    df['PLEsoil_day'] = (df.s_day*df.Asoil_day+df.rho_day*Cp*(1.0-df.Fc)*df.vpd_day*df.daylength/df.ras) * \
                        (1.0-df.Fwet_day)/(df.s_day+df.gama_day*df.rtot_day/df.ras)

    return df


# total soil LE
def LEsoil_day(df, sm_open, sm_close):

    df = PLE_soil_day(df)
    df = LEwetsoil_day(df)

    if sm_open is not None and sm_close is not None:
        df['fSM_day'] = df.fSM
    else:
        df['fSM_day'] = (df.rh*0.01) ** (df.vpd_day/250.0)

    df['LEsoil_day'] = (df.PLEsoil_day*df.fSM_day)+df.LEwetsoil_day
    return df


# total soil LE
def LEsoil_night(df, sm_open, sm_close):

    df = PLE_soil_night(df)
    df = LEwetsoil_night(df)

    if sm_open is not None and sm_close is not None:
        df['fSM_night'] = df.fSM
    else:
        df['fSM_night'] = (df.rh*0.01) ** (df.vpd_night/250.0)

    df['LEsoil_night'] = (df.PLEsoil_night*df.fSM_night)+df.LEwetsoil_night
    return df


# -------------------------------------------------Evapotranspiration------------------------------------------
# calculate total ET
def MOD16(df, tmin_open, tmin_close, vpd_close, vpd_open, gl_sh,
          gl_e_wv, Cl, rbl_min, rbl_max, sm_open=None, sm_close=None):

    df = LEwetfun_day(df, gl_sh=gl_sh, gl_e_wv=gl_e_wv)
    df = LEwetfun_night(df, gl_sh=gl_sh, gl_e_wv=gl_e_wv)

    df = LEtrans_day(df, Cl=Cl, gl_sh=gl_sh, vpd_open=vpd_open, vpd_close=vpd_close,
                     tmin_open=tmin_open, tmin_close=tmin_close, sm_open=sm_open, sm_close=sm_close)
    df = LEtrans_night(df, gl_sh=gl_sh)
    df = PLE_plant_day(df)
    df = PLE_plant_night(df)

    df = rtot_day(df, rbl_max=rbl_max, rbl_min=rbl_min, vpd_close=vpd_close, vpd_open=vpd_open)
    df = rtot_night(df, rbl_max=rbl_max, rbl_min=rbl_min, vpd_close=vpd_close, vpd_open=vpd_open)

    df = LEwetsoil_day(df)
    df = LEwetsoil_night(df)

    df = PLE_soil_day(df)
    df = PLE_soil_night(df)

    df = LEsoil_day(df, sm_open=sm_open, sm_close=sm_close)
    df = LEsoil_night(df, sm_open=sm_open, sm_close=sm_close)

    df = df.fillna(0)

    df['LE_day'] = df.LEwet_day + df.LEtrans_day + df.LEsoil_day
    df['LE_night'] = df.LEwet_night + df.LEtrans_night + df.LEsoil_night

    df['ET_day'] = df.LE_day/df.lamda_day
    df['ET_night'] = df.LE_night/df.lamda_night
    df['LE'] = df.LE_day + df.LE_night
    df['simulation'] = df.ET_day + df.ET_night

    return df
