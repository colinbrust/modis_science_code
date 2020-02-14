import ee
from msc.utils import Bplut as bp, GetModisFpar as fp


# Global constants
time = 'system:time_start'
day = (24 * 60 * 60 * 1000)


def MOD16(roi: ee.Geometry, year: int, **kwargs) -> ee.ImageCollection:
    """
    :param roi: ee.Geometry of region to run the model across
    :param year: Year to run the model in
    :param kwargs: Additional arguments including:
        meteorology: ee.ImageCollection, meteorology containing bands named 'vpd', 'srad', 'tmmn', 'tmmx',
            'rmin', 'rmax', 'vpd'. ***Required***
        daylength: ee.ImageCollection, gridded product giving daylength in seconds ***Required***
        elev: ee.Image, elevation in meters ***Required***
        LAI: ee.ImageCollection of LAI,
        Fc: ee.ImageCollection of FPAR/Fc,
        smap_sm: ee.ImageCollection with 'rzMean' and 'fSM' bands returned from msc.utils.Smap.get_smap.
        bplut: ee.Image of bplut params returned from a function in msc.utils.Bplut
        all_variables: bool, whether to export just ET or all intermediate variables as well.
    :return: ee.ImageCollection containing modeled evapotranspiration
    """
    try:
        meteorology = ee.ImageCollection(kwargs.pop('meteorology'))
        daylength = kwargs.pop('daylength')
        elev = kwargs.pop('elev')
    except KeyError as e:
        print('Missing required {} input. Rerun model including this argument.'.format(e))

    start = ee.Date.fromYMD(year, 1, 1)
    end = ee.Date.fromYMD(year + 1, 1, 1)

    # Function to join the daily daylength, albedo and meteorology data together
    def dataJoin(left, right):
        filt = ee.Filter.equals(
            leftField=time,
            rightField=time)

        join = ee.Join.saveBest(
            matchKey='match',
            measureKey='delta_t')

        return ee.ImageCollection(join.apply(left, right, filt)) \
            .map(lambda img: img.addBands(img.get('match')))

    if 'LAI' in kwargs and 'Fc' in kwargs:
        fp_lai = dataJoin(kwargs.pop('LAI'), kwargs.pop('Fc'))
        fp_lai = fp.modis_fpar_lai(roi, year, fp_lai)
    elif ('LAI' in kwargs and 'Fc' not in kwargs) or ('LAI' not in kwargs and 'Fc' in kwargs):
        raise ValueError('If either LAI or Fc is used as a kwarg, the other must also be included as an argument.\n'
                         'Please include both LAI and Fc.')
    else:
        fp_lai = fp.modis_fpar_lai(roi, year)

    LAI = ee.ImageCollection(fp_lai.get('LAI'))
    Fc = ee.ImageCollection(fp_lai.get('FPAR'))

    proj = ee.Image(LAI.first()).projection()

    # Function to downscale inputs to match MODIS projection and resolution
    def match_proj(img):
        img = img.resample('bilinear').reproject(
            crs=proj.crs(),
            scale=500
        ).copyProperties(img, ['system:time_start', 'system:index'])

        return img

    # Import code that contains a spatial BPLUT
    bplut = bp.m16_bplut(roi, start, end) if 'bplut' not in kwargs else kwargs.pop('bplut')

    # Filter meteorological data
    meteorology = meteorology \
        .filterDate(start, end) \
        .map(lambda img: img.clip(roi)) \
        .map(match_proj)

    def avg_temp(img):
        tavg = img.expression('(tmin+tmax)/2-273.15', {
            'tmin': img.select('tmmn'),
            'tmax': img.select('tmmx')
        }).rename('Tavg')

        return img.addBands(tavg)

    meteorology = meteorology.map(lambda img: avg_temp(img))

    # import daylength data
    daylength = daylength \
        .select('dayl') \
        .filterDate(start, end) \
        .map(lambda img: img.clip(roi))

    # QA/QC filter for MODIS albedo
    def albedoqc(img):
        qc = img.select(['BRDF_Albedo_Band_Mandatory_Quality_shortwave']) \
            .bitwiseAnd(1).eq(0)

        return img.updateMask(qc)

    # Import MODIS albedo data, daily, 500m
    albedo = ee.ImageCollection('MODIS/006/MCD43A3') \
        .filterBounds(roi) \
        .filterDate(start.advance(-1, 'month'), end) \
        .map(albedoqc)

    albedo = albedo.select('Albedo_WSA_shortwave') \
        .map(lambda img: img.rename('albedo').clip(roi))

    # interpolate albedo with using the mean annual albedo each year.
    mean_albedo = albedo.mean()

    def albedo_fun(img):
        masked = img.unmask()
        filled = masked.where(masked.eq(0), mean_albedo)

        return filled.rename('albedo') \
            .copyProperties(img, ['system:time_start']) \
            .copyProperties(img, ['system:time_end'])

    albedo_interp = albedo.map(albedo_fun)


    def dataJoin2(left, right):
        filt = ee.Filter.maxDifference(
            difference=day,
            leftField=time,
            rightField=time)

        join = ee.Join.saveBest(
            matchKey='match',
            measureKey='delta_t')

        return ee.ImageCollection(join.apply(left, right, filt)) \
            .map(lambda img: img.addBands(img.get('match')))

    # Group all input data into one imageCollection
    meteo = dataJoin(daylength, albedo_interp)
    meteo = dataJoin2(meteo, meteorology)
    meteo = dataJoin2(meteo, Fc)
    meteo = dataJoin2(meteo, LAI)

    if 'smap_sm' in kwargs:
        smap_sm = kwargs['smap_sm']
        smap_sm = smap_sm.map(match_proj)
        meteo = dataJoin2(meteo, smap_sm)

    meteo = meteo.map(lambda img: img.clip(roi))
    meteo = meteo.filterDate(start, end)

    # calculate annual mean daily temperature in degree C and add to bplut for soil heat flux calculations
    Tannual = meteo.select('Tavg').mean().rename('Tannual')
    bplut = bplut.addBands(Tannual)

    # Constants necessary for ET calculation
    SBC = 5.67 / 1e8  # (W/(m2 K4)) Stefan-Boltzmann constant
    Cp = 1013.0
    epsl = 0.622

    # Function to calculate net radiation
    def calc_rn(temp, albedo, swrad, daylength):
        t = temp.expression('pow((273.15+Tavg),4)', {  # temp: degree C
            'Tavg': temp})

        ea = temp.expression('1-0.26*exp((-7.77)*0.0001*T*T)', {
            'T': temp})

        Rnet = t.expression('((1-a)*R+(ea-0.97)*t*sbc)*d', {  # Rn unit:J/m2
            'a': albedo,
            'R': swrad,  # W/m2
            't': t,
            'ea': ea,
            'd': daylength,
            'sbc': SBC
        })
        return Rnet  # J/day/m2W

    # Function to calculate soil heat flux
    def calc_soil_heat_flux(bplut, Rn, temp, daylength, tday, tnight):
        gsoil = temp.expression('(4.73*T-20.87)*day', {
            'T': temp,
            'day': daylength
        })

        con1x = bplut.select('Tannual').gte(bplut.select('tmin_close')) \
            .And(bplut.select('Tannual').lt(ee.Image(25))) \
            .And((tday.subtract(tnight)).gte(ee.Image(5)))
        con1 = ee.Image(0)
        con1 = con1.where(con1x, 1)

        gsoil = gsoil.multiply(con1)
        con2 = gsoil.expression('(abs(g) > (0.39*abs(a))) ? 1 : 0', {
            'g': gsoil,
            'a': Rn
        })

        g2 = Rn.multiply(0.39)
        gsoil = gsoil.where(con2, g2)
        return gsoil

    # Functions to caculate energy available at canopy and soil surfaces
    def calc_canopy_rad(Rn, Fc):
        Ac = Rn.expression('A*Fc', {
            'A': Rn,
            'Fc': Fc
        })
        return Ac

    def calc_soil_rad(Rn, Fc, G):
        Asoil = Rn.expression('(1-Fc)*A-G', {
            'Fc': Fc,
            'A': Rn,
            'G': G
        })
        return Asoil

    # Function to calculate fraction of soil surface that is saturated
    def calc_wet_soil(rh):  # rh unit %
        value = rh.expression('pow((rh*0.01),4)', {'rh': rh})
        Fwet = ee.Image(0)
        cond1 = rh.gte(ee.Image(70)).And(rh.lte(ee.Image(100)))
        Fwet = Fwet.where(cond1, value)
        return Fwet.rename('Fwet')

        # Function to calculate atmospheric pressure given elevation

    def Pafun(elevation):
        t1 = elevation.expression('1.0-(LRstd*elev)/Tstd', {
            'LRstd': 0.0065,
            'elev': elevation,
            'Tstd': 288.15
        })
        t2 = elevation.expression('Gstd/(LRstd*RR/MA)', {
            'Gstd': 9.80665,
            'LRstd': 0.0065,
            'RR': 8.3143,
            'MA': 0.0289644
        })
        Pa = t1.expression('Pstd*pow(t1,t2)', {
            'Pstd': 101325.0,
            't1': t1,
            't2': t2
        })
        return Pa.reproject(proj)

    # 1 pa = 0.01 mbar
    # unit of rho:kg/m3, T: degree C, rh:%
    # Air density calculation
    # air density calculation:Buoyancy Correction and Air Density Measurement. National Physical Laboratory
    # http:#www.npl.co.uk/upload/pdf/buoycornote.pdf
    def calc_rho(rh, pa, temp):
        rho = rh.expression('(0.34844*P*0.01-rh*(0.00252*T-0.020582))/(T+273.15)', {
            'P': pa,
            'rh': rh,
            'T': temp
        })
        return rho

    # Function to calculate saturated vapor pressure
    def calc_saturated_vapor_pressure(temp):
        svp = temp.expression('610.7*exp(17.38*t/(t+239))', {
            't': temp
        })
        return svp  # Pa

    # Function to calculate latent heat of vaporization
    def calc_latent_heat_vaporization(temp):
        lamda = temp.expression('(2.501-0.002361*t)*1e6', {
            't': temp
        })
        return lamda

    # Function to calculate slope of saturation water vapor pressure curve
    def calc_svp_slope(svp, temp):  # unit: Pa/k
        s = svp.expression('17.38*239.0*svp/(t*t)', {
            'svp': svp,
            't': temp.add(239.0)
        })
        return s

    # Function to calculate latent heat flux for wet canopy
    def calc_wet_canopy_evap(Ac, Fc, rho, s, vpd, daylength, Fwet, Pa, lamda, bplut, LAI, temp):
        rhc = LAI.expression('1.0/(gl*lai*fwet)', {
            'gl': bplut.select('gl_sh'),
            'lai': LAI,
            'fwet': Fwet
        })
        rrc = temp.expression('rho*Cp/(4.0*SBC*t*t*t)', {
            'rho': rho,
            'Cp': Cp,
            'SBC': SBC,
            't': temp.add(273.15)
        })
        rhrc = rhc.expression('rhc*rrc/(rhc+rrc)', {
            'rhc': rhc,
            'rrc': rrc
        })
        rvc = LAI.expression('1.0/(gl_e_wv*lai*fwet)', {
            'gl_e_wv': bplut.select('gl_e_wv'),
            'lai': LAI,
            'fwet': Fwet
        })

        LEwet_c = Ac.expression('(s*Ac+rho*Cp*vpd*Fc*daylength/rhrc)*Fwet/(s+(Pa*Cp*rvc/(lamda*epsl*rhrc)))', {
            's': s,
            'Ac': Ac,
            'Fc': Fc,
            'rho': rho,
            'Cp': Cp,
            'vpd': vpd,
            'daylength': daylength,
            'rhrc': rhrc,
            'Fwet': Fwet,
            'Pa': Pa,
            'rvc': rvc,
            'lamda': lamda,
            'epsl': epsl
        })
        return LEwet_c

    # Function to calculate PFT dependent VPD constraint scalar
    def vpd_scalar(bplut, vpd):
        mvpd = vpd.expression('(vpdclose-vpd)/(vpdclose-vpdopen)', {
            'vpdclose': bplut.select('vpd_close'),
            'vpd': vpd,
            'vpdopen': bplut.select('vpd_open')
        }).clamp(0, 1)

        return mvpd

    # Function to calculate PFT dependent Tmin constraint scalar
    def t_scalar(bplut, temp):
        mT = temp.expression('(t-tminclose)/(tminopen-tminclose)', {
            't': temp,
            'tminclose': bplut.select('tmin_close'),
            'tminopen': bplut.select('tmin_open')
        }).clamp(0, 1)
        return mT

    def sm_scalar(bplut, sm):
        mSM = sm.expression('(sm-smclose)/(smopen-smclose)', {
            'sm': sm,
            'smclose': bplut.select('sm_close'),
            'smopen': bplut.select('sm_open')}).clamp(0, 1)
        return mSM

    # Function to calculate psychrometric constant
    def calc_psychrometric_const(pa, lamda):
        gama = pa.expression('Cp*pa/(lamda*epsl)', {
            'Cp': Cp,
            'pa': pa,
            'lamda': lamda,
            'epsl': epsl
        })
        return gama

    # Function to calculate transpiration that is used in both daytime and nighttime calculations.
    def le_trans_base(rcorr, Gs1, temp, bplut, LAI, Fwet, rho, gama, s, Ac, Fc, vpd, daylength):

        Gcu = rcorr.multiply(0.00001)
        Gs2 = bplut.select('gl_sh')
        Cc = LAI.expression('lai*(1.0-fwet)*(Gs2*(Gs1+Gcu))/(Gs1+Gs2+Gcu)', {
            'lai': LAI,
            'fwet': Fwet,
            'Gs1': Gs1,
            'Gs2': Gs2,
            'Gcu': Gcu
        })

        cond1 = LAI.gt(ee.Image(0)).And((ee.Image(1).subtract(Fwet)).gt(ee.Image(0)))
        cond = ee.Image(0)
        cond = cond.where(cond1, 1)
        Cc = Cc.multiply(cond)

        rs = Cc.expression('1.0/Cc', {'Cc': Cc})

        # calculate aerodynamic resistance
        rh = bplut.expression('1.0/gl_sh', {'gl_sh': bplut.select('gl_sh')})
        rr = rho.expression('rho*Cp/(4.0*SBC*t*t*t)', {
            'rho': rho,
            'Cp': Cp,
            'SBC': SBC,
            't': temp.add(273.15)
        })
        ra = rr.expression('(rh*rr)/(rh+rr)', {'rh': rh, 'rr': rr})

        LEtrans = Ac.expression('(s*Ac+rho*Cp*vpd*Fc*daylength/ra)*(1-fwet)/(s+gama*(1+rs/ra))', {
            's': s,
            'Ac': Ac,
            'Fc': Fc,
            'rho': rho,
            'Cp': Cp,
            'vpd': vpd,
            'daylength': daylength,
            'ra': ra,
            'fwet': Fwet,
            'gama': gama,
            'rs': rs
        })
        return LEtrans

    # Function to calculate daytime plant transpiration
    def calc_le_trans_day(temp, pa, bplut, LAI, Fwet, rho, gama, s, Ac, Fc, vpd, daylength, sm=None):
        ta = temp.expression('pow(((t+273.15)/293.15),1.75)', {'t': temp})
        rcorr = pa.expression('1.0/((101300/pa)*ta)', {
            'pa': pa,
            'ta': ta
        })
        mvpd = vpd_scalar(bplut, vpd)
        mT = t_scalar(bplut, temp)
        if 'smap_sm' in kwargs:
            mSM = sm_scalar(bplut, sm)
            Gs1 = bplut.expression('Cl*mvpd*mT*mSM*rcorr', {
                'Cl': bplut.select('Cl'),
                'mvpd': mvpd,
                'mT': mT,
                'mSM': mSM,
                'rcorr': rcorr
            })
        else:
            Gs1 = bplut.expression('Cl*mvpd*mT*rcorr', {
                'Cl': bplut.select('Cl'),
                'mvpd': mvpd,
                'mT': mT,
                'rcorr': rcorr
            })

        return le_trans_base(rcorr=rcorr, Gs1=Gs1, temp=temp, bplut=bplut, LAI=LAI, Fwet=Fwet, rho=rho, gama=gama,
                             s=s, Ac=Ac, Fc=Fc, vpd=vpd, daylength=daylength)

    # Function to calculate night time plant tranpiration
    def calc_le_trans_night(temp, pa, bplut, LAI, Fwet, rho, gama, s, Ac, Fc, vpd, daylength):
        ta = temp.expression('pow(((t+273.15)/293.15),1.75)', {'t': temp})
        rcorr = pa.expression('1.0/((101300/pa)*ta)', {
            'pa': pa,
            'ta': ta
        })

        Gs1 = ee.Image(0)

        return le_trans_base(rcorr=rcorr, Gs1=Gs1, temp=temp, bplut=bplut, LAI=LAI, Fwet=Fwet, rho=rho, gama=gama,
                             s=s, Ac=Ac, Fc=Fc, vpd=vpd, daylength=daylength)

    # Function to calculate total soil resistance to water vapor transport
    def calc_total_soil_resistance(temp, vpd, bplut, pa):
        ta = temp.expression('pow(((t+273.15)/293.15),1.75)', {'t': temp})
        rcorr = ta.expression('1.0/((101300.0/pa)*ta)', {
            'pa': pa,
            'ta': ta
        })

        cond1 = ee.Image(0)
        cond2 = ee.Image(0)
        cond3 = ee.Image(0)
        cond1x = vpd.lte(bplut.select('vpd_open'))
        cond1 = cond1.where(cond1x, 1)
        cond2x = vpd.gte(bplut.select('vpd_close'))
        cond2 = cond2.where(cond2x, 1)
        cond3x = vpd.gt(bplut.select('vpd_open')).And(vpd.lt(bplut.select('vpd_close')))
        cond3 = cond3.where(cond3x, 1)

        rt1 = bplut.select('rbl_max').multiply(cond1)
        rt2 = bplut.select('rbl_min').multiply(cond2)
        rt3 = bplut.expression('rblmax-((rblmax-rblmin)*(vpdclose-vpd))/(vpdclose-vpdopen)', {
            'rblmax': bplut.select('rbl_max'),
            'rblmin': bplut.select('rbl_min'),
            'vpdclose': bplut.select('vpd_close'),
            'vpd': vpd,
            'vpdopen': bplut.select('vpd_open')
        })

        rt3 = rt3.multiply(cond3)

        rtotc = rt1.add(rt2).add(rt3)
        rtot = rcorr.multiply(rtotc)
        return rtot

    # Function to calculate resistance to radiative heat transfer
    def calc_rht_resistance(rho, temp):
        rrs = rho.expression('rho*Cp/(4.0*SBC*T*T*T)', {
            'rho': rho,
            'Cp': Cp,
            'SBC': SBC,
            'T': temp.add(273.15)
        })
        return rrs

    # Function to calculate latent heat flux at saturated soil surface
    def calc_sat_soil_evap(rtot, rrs, rho, vpd, s, Asoil, Fc, Fwet, gama, daylength):
        ras = rrs.expression('rhs*rrs/(rhs+rrs)', {
            'rhs': rtot,
            'rrs': rrs
        })
        LEwet_soil = s.expression('(s*Asoil+rho*Cp*(1.0-Fc)*vpd*daylength/ras)*fwet/(s+gama*rtot/ras)', {
            's': s,
            'Asoil': Asoil,
            'rho': rho,
            'Cp': Cp,
            'Fc': Fc,
            'vpd': vpd,
            'daylength': daylength,
            'ras': ras,
            'fwet': Fwet,
            'gama': gama,
            'rtot': rtot,
        })
        return LEwet_soil

    # Function to calculate potential soil evaporation
    def calc_potential_soil_evap(rtot, rrs, s, Asoil, rho, Fc, vpd, Fwet, daylength, gama):
        ras = rrs.expression('rhs*rrs/(rhs+rrs)', {
            'rhs': rtot,
            'rrs': rrs
        })

        PLEsoil = s.expression('(s*Asoil+rho*Cp*(1.0-Fc)*vpd*daylength/ras)*(1.0-fwet)/(s+gama*rtot/ras)', {
            's': s,
            'Asoil': Asoil,
            'rho': rho,
            'Cp': Cp,
            'Fc': Fc,
            'vpd': vpd,
            'daylength': daylength,
            'ras': ras,
            'fwet': Fwet,
            'gama': gama,
            'rtot': rtot
        })
        return PLEsoil

    # Function to calculate latent heat flux for non-saturated soil
    def calc_total_soil_evap(rh, vpd, LEwetsoil, PLEsoil, rew):
        
        if 'smap_sm' in kwargs:
            return PLEsoil.multiply(rew).add(LEwetsoil)
        else:
            fSM = rh.expression('pow((rh*0.01),(vpd/beta))', {
                'rh': rh,
                'vpd': vpd,
                'beta': 250
            })
            return PLEsoil.multiply(fSM).add(LEwetsoil)

    # Aggregate all functions to calculate total ET
    def calc_et(img):
        swrad = img.select('srad')            # W/m2
        albedo = img.select('albedo').multiply(0.001)
        ta = img.select('Tavg')
        ta_day = ta   # degree C
        ta_night = img.select('tmmn').subtract(273.15)
        rhmax = img.select('rmax')
        rhmin = img.select('rmin')               # rh_night
        rh = (rhmax.add(rhmin)).multiply(0.5)   # % rh_day
        vpd_day = img.select('vpd').multiply(1000)                 # Pa
        daylength = img.select('dayl')
        nightlength = daylength.expression('24*3600.0-dayl', {'dayl': daylength})
        Fc = img.select('Fc')
        LAI = img.select('LAI')
        svp_night = calc_saturated_vapor_pressure(ta_night)   # Pa
        svp_day = calc_saturated_vapor_pressure(ta_day)
        vpd_night = svp_night.multiply(ee.Image(1).subtract((rhmin.multiply(0.01))))  # Pa
        pa = Pafun(elev)

        doy = ee.Date(img.get('system:time_start'))

        sm_rz, rew = (img.select('rzMean'), img.select('fSM')) if 'smap_sm' in kwargs else (ee.Image(0), ee.Image(0))

        # calculate each parameter values in day and night, respectively
        Rn_day = calc_rn(temp=ta_day, albedo=albedo, swrad=swrad, daylength=daylength)  # J/day/m2
        Rn_night = calc_rn(temp=ta_night, albedo=albedo, swrad=0, daylength=nightlength)

        Gsoil_day = calc_soil_heat_flux(bplut=bplut, Rn=Rn_day, temp=ta_day, daylength=daylength, tday=ta_day,
                                        tnight=ta_night)
        Gsoil_night = calc_soil_heat_flux(bplut=bplut, Rn=Rn_night, temp=ta_night, daylength=nightlength, tday=ta_day,
                                          tnight=ta_night)

        #  following the User Guide (2017 version)
        G_day = Gsoil_day.multiply(ee.Image(1).subtract(Fc))
        cond1 = Rn_day.lt(G_day)
        G_day = G_day.where(cond1, Rn_day)

        G_night = Gsoil_night.multiply(ee.Image(1).subtract(Fc))
        cond2 = (Rn_night.subtract(G_night)).lt(Rn_day.multiply(-0.5))
        g1 = Rn_night.add(Rn_day.multiply(0.5))
        G_night = G_night.where(cond2, g1)

        Ac_day = calc_canopy_rad(Rn=Rn_day, Fc=Fc)
        Ac_night = calc_canopy_rad(Rn=Rn_night, Fc=Fc)
        Asoil_day = calc_soil_rad(Rn=Rn_day, Fc=Fc, G=G_day)
        Asoil_night = calc_soil_rad(Rn=Rn_night, Fc=Fc, G=G_night)

        Fwet_day = calc_wet_soil(rh=rh)
        Fwet_night = calc_wet_soil(rh=rhmin)

        rho_day = calc_rho(rh=rh, pa=pa, temp=ta_day)
        rho_night = calc_rho(rh=rhmin, pa=pa, temp=ta_night)

        lamda_day = calc_latent_heat_vaporization(temp=ta_day)
        lamda_night = calc_latent_heat_vaporization(temp=ta_night)

        s_day = calc_svp_slope(svp=svp_day, temp=ta_day)
        s_night = calc_svp_slope(svp=svp_night, temp=ta_night)

        LEwet_day = calc_wet_canopy_evap(Ac=Ac_day, Fc=Fc, rho=rho_day, s=s_day, vpd=vpd_day, daylength=daylength,
                                         Fwet=Fwet_day, Pa=pa, lamda=lamda_day, bplut=bplut, LAI=LAI, temp=ta_day)
        LEwet_night = calc_wet_canopy_evap(Ac=Ac_night, Fc=Fc, rho=rho_night, s=s_night, vpd=vpd_night,
                                           daylength=nightlength, Fwet=Fwet_night, Pa=pa, lamda=lamda_night,
                                           bplut=bplut, LAI=LAI, temp=ta_night)

        gama_day = calc_psychrometric_const(pa=pa, lamda=lamda_day)
        gama_night = calc_psychrometric_const(pa=pa, lamda=lamda_night)

        LEtrans_day = calc_le_trans_day(temp=ta_day, pa=pa, bplut=bplut, LAI=LAI, Fwet=Fwet_day, rho=rho_day,
                                        gama=gama_day, s=s_day, Ac=Ac_day, Fc=Fc, vpd=vpd_day, daylength=daylength,
                                        sm=sm_rz)
        LEtrans_night = calc_le_trans_night(temp=ta_night, pa=pa, bplut=bplut, LAI=LAI, Fwet=Fwet_night, rho=rho_night,
                                            gama=gama_night, s=s_night, Ac=Ac_night, Fc=Fc, vpd=vpd_night,
                                            daylength=nightlength)

        rtot_day = calc_total_soil_resistance(temp=ta_day, vpd=vpd_day, bplut=bplut, pa=pa)
        rtot_night = calc_total_soil_resistance(temp=ta_night, vpd=vpd_night, bplut=bplut, pa=pa)
        rrs_day = calc_rht_resistance(rho=rho_day, temp=ta_day)
        rrs_night = calc_rht_resistance(rho=rho_night, temp=ta_night)

        LEwetsoil_day = calc_sat_soil_evap(rtot=rtot_day, rrs=rrs_day, rho=rho_day, vpd=vpd_day, s=s_day,
                                           Asoil=Asoil_day, Fc=Fc, Fwet=Fwet_day, gama=gama_day, daylength=daylength)
        LEwetsoil_night = calc_sat_soil_evap(rtot=rtot_night, rrs=rrs_night, rho=rho_night, vpd=vpd_night, s=s_night,
                                             Asoil=Asoil_night, Fc=Fc, Fwet=Fwet_night, gama=gama_night,
                                             daylength=nightlength)
        PLEsoil_day = calc_potential_soil_evap(rtot=rtot_day, rrs=rrs_day, s=s_day, Asoil=Asoil_day, rho=rho_day, Fc=Fc,
                                               vpd=vpd_day, Fwet=Fwet_day, daylength=daylength, gama=gama_day)
        PLEsoil_night = calc_potential_soil_evap(rtot=rtot_night, rrs=rrs_night, s=s_night, Asoil=Asoil_night,
                                                 rho=rho_night, Fc=Fc, vpd=vpd_night, Fwet=Fwet_night,
                                                 daylength=nightlength, gama=gama_night)
        LEsoil_day = calc_total_soil_evap(rh=rh, vpd=vpd_day, LEwetsoil=LEwetsoil_day, PLEsoil=PLEsoil_day, rew=rew)
        LEsoil_night = calc_total_soil_evap(rh=rh, vpd=vpd_night, LEwetsoil=LEwetsoil_night, PLEsoil=PLEsoil_night,
                                            rew=rew)

        # calculate total daily LE, PLE, ET and PET.
        LE_day = LEwet_day.add(LEtrans_day).add(LEsoil_day)
        ET_day = LE_day.divide(lamda_day)
        LE_night = LEwet_night.add(LEtrans_night).add(LEsoil_night)
        ET_night = LE_night.divide(lamda_night)
        LE = LE_day.add(LE_night)
        ET = ET_day.add(ET_night)

        img_out = Rn_day.addBands(Rn_night).addBands(G_day).addBands(G_night).addBands(Ac_day).addBands(Ac_night)\
            .addBands(Asoil_day).addBands(Asoil_night).addBands(Fwet_day).addBands(Fwet_night).addBands(rho_day)\
            .addBands(rho_night).addBands(lamda_day).addBands(lamda_night).addBands(s_day).addBands(s_night)\
            .addBands(LEwet_day).addBands(LEwet_night).addBands(gama_day).addBands(gama_night).addBands(LEtrans_day)\
            .addBands(LEtrans_night).addBands(rtot_day).addBands(rtot_night).addBands(rrs_day).addBands(rrs_night)\
            .addBands(LEwetsoil_day).addBands(LEwetsoil_night).addBands(PLEsoil_day).addBands(PLEsoil_night)\
            .addBands(LEsoil_day).addBands(LEsoil_night).addBands(LE).addBands(ET).addBands(daylength) .addBands(swrad)\
            .addBands(albedo).addBands(ta).addBands(ta_day).addBands(ta_night).addBands(rhmax).addBands(rhmin)\
            .addBands(rh).addBands(vpd_day).addBands(nightlength).addBands(Fc).addBands(LAI).addBands(svp_night) \
            .addBands(svp_day).addBands(vpd_night).addBands(pa)\
            .select([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                     18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36,
                     37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50],
                    ['Rn_day', 'Rn_night', 'G_day', 'G_night', 'Ac_day', 'Ac_night', 'Asoil_day', 'Asoil_night',
                     'Fwet_day', 'Fwet_night', 'rho_day', 'rho_night', 'lamda_day', 'lamda_night', 's_day',
                     's_night', 'LEwet_day', 'LEwet_night', 'gama_day', 'gama_night', 'LEtrans_day', 'LEtrans_night',
                     'rtot_day', 'rtot_night', 'rrs_day', 'rrs_night', 'LEwetsoil_day', 'LEwetsoil_night',
                     'PLEsoil_day', 'PLEsoil_night', 'LEsoil_day', 'LEsoil_night', 'LE', 'ET', 'daylength',
                     'swrad', 'albedo', 'ta', 'ta_day', 'ta_night', 'rhmax', 'rhmin', 'rh', 'vpd_day',
                     'nightlength', 'Fc', 'LAI', 'svp_night', 'svp_day', 'vpd_night', 'pa']) \
            .copyProperties(img, [time]) \
            .set({'Date': doy})

        if 'smap_sm' in kwargs:
            img_out = ee.Image(img_out).addBands(sm_rz).addBands(rew).addBands(img.select('surfMean'))

        return ee.Image(img_out)

    et_out = meteo.map(calc_et)

    all_variables = kwargs.pop('all_variables') if 'all_variables' in kwargs else False
    if not all_variables:
        return et_out.select('ET')
    return et_out
