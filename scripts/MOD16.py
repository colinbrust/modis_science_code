import ee
from scripts import GetModisFpar as fp
from scripts import BPLUT as bp

ee.Initialize()

# Global constants
time = 'system:time_start'
day = (24 * 60 * 60 * 1000)

def MOD16(roi, year, meteorology, daylength, elev):
    start = ee.Date.fromYMD(year, 1, 1)
    end = ee.Date.fromYMD(year + 1, 1, 1)

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
    bplut = bp.m16_bplut(roi, start, end)

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

    # Function to join the daily daylength, albedo and meteorology data together
    def dataJoin(left, right):
        filter = ee.Filter.equals(
            leftField=time,
            rightField=time)

        join = ee.Join.saveBest(
            matchKey='match',
            measureKey='delta_t')

        return ee.ImageCollection(join.apply(left, right, filter)) \
            .map(lambda img: img.addBands(img.get('match')))

    def dataJoin2(left, right):
        filter = ee.Filter.maxDifference(
            difference=day,
            leftField=time,
            rightField=time)

        join = ee.Join.saveBest(
            matchKey='match',
            measureKey='delta_t')

        return ee.ImageCollection(join.apply(left, right, filter)) \
            .map(lambda img: img.addBands(img.get('match')))

    # Group all input data into one imageCollection
    meteo = dataJoin(daylength, albedo_interp)
    meteo = dataJoin2(meteo, meteorology)
    meteo = dataJoin2(meteo, Fc)
    meteo = dataJoin2(meteo, LAI)
    meteo = meteo.map(lambda img: img.clip(roi))

    # calculate annual mean daily temperature in degree C and add to bplut for soil heat flux calculations
    Tannual = meteo.select('Tavg').mean().rename('Tannual')
    bplut = bplut.addBands(Tannual)

    # Constants necessary for ET calculation
    SBC = 5.67 / 1e8  # (W/(m2 K4)) Stefan-Boltzmann constant
    Cp = 1013.0
    epsl = 0.622

    # Function to calculate net radiation
    def Rnfun(temp, albedo, swrad, daylength):
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
    def Gsoilfun(bplut, Rn, temp, daylength, tday, tnight, fc):
        gsoil = temp.expression('(4.73*T-20.87)*day', {
            'T': temp,
            'day': daylength
        })

        con1x = bplut.select('Tannual').gte(bplut.select('Tminclose')) \
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
    def Acfun(Rn, Fc):
        Ac = Rn.expression('A*Fc', {
            'A': Rn,
            'Fc': Fc
        })
        return Ac

    def Asoilfun(Rn, Fc, G):
        Asoil = Rn.expression('(1-Fc)*A-G', {
            'Fc': Fc,
            'A': Rn,
            'G': G
        })
        return Asoil

    # Function to calculate fraction of soil surface that is saturated
    def Fwetfun(rh):  # rh unit %
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
    def rhofun(rh, pa, temp):
        rho = rh.expression('(0.34844*P*0.01-rh*(0.00252*T-0.020582))/(T+273.15)', {
            'P': pa,
            'rh': rh,
            'T': temp
        })
        return rho

    # Function to calculate saturated vapor pressure
    def svpfun(temp):
        svp = temp.expression('610.7*exp(17.38*t/(t+239))', {
            't': temp
        })
        return svp  # Pa

    # Function to calculate latent heat of vaporization
    def lamdafun(temp):
        lamda = temp.expression('(2.501-0.002361*t)*1e6', {
            't': temp
        })
        return lamda

    # Function to calculate slope of saturation water vapor pressure curve
    def sfun(svp, temp):  # unit: Pa/k
        s = svp.expression('17.38*239.0*svp/(t*t)', {
            'svp': svp,
            't': temp.add(239.0)
        })
        return s

    # Function to calculate latent heat flux for wet canopy
    def LEwetfun(Ac, Fc, rho, s, vpd, daylength, Fwet, Pa, lamda, bplut, LAI, temp):
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
            'vpdclose': bplut.select('VPDclose'),
            'vpd': vpd,
            'vpdopen': bplut.select('VPDopen')
        }).clamp(0, 1)

        return mvpd

    # Function to calculate PFT dependent Tmin constraint scalar
    def T_scalar(bplut, temp):
        mT = temp.expression('(t-tminclose)/(tminopen-tminclose)', {
            't': temp,
            'tminclose': bplut.select('Tminclose'),
            'tminopen': bplut.select('Tminopen')
        }).clamp(0, 1)
        return mT

    # Function to calculate psychrometric constant
    def gamafun(pa, lamda):
        gama = pa.expression('Cp*pa/(lamda*epsl)', {
            'Cp': Cp,
            'pa': pa,
            'lamda': lamda,
            'epsl': epsl
        })
        return gama

    # Function to calculate daytime plant transpiration
    def LEtransfun(temp, pa, bplut, LAI, Fwet, rho, lamda, gama, s, Ac, Fc, vpd, daylength):
        ta = temp.expression('pow(((t+273.15)/293.15),1.75)', {'t': temp})
        rcorr = pa.expression('1.0/((101300/pa)*ta)', {
            'pa': pa,
            'ta': ta
        })
        mvpd = vpd_scalar(bplut, vpd)
        mT = T_scalar(bplut, temp)
        Gs1 = bplut.expression('Cl*mvpd*mT*rcorr', {
            'Cl': bplut.select('Cl'),
            'mvpd': mvpd,
            'mT': mT,
            'rcorr': rcorr
        })
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

    # Function to calculate night time plant tranpiration
    def LEtransfun2(temp, pa, bplut, LAI, Fwet, rho, lamda, gama, s, Ac, Fc, vpd, daylength):
        ta = temp.expression('pow(((t+273.15)/293.15),1.75)', {'t': temp})
        rcorr = pa.expression('1.0/((101300/pa)*ta)', {
            'pa': pa,
            'ta': ta
        })
        mvpd = vpd_scalar(bplut, vpd)
        mT = T_scalar(bplut, temp)
        Gs1 = ee.Image(0)
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

    # Potential plant transpiration calculated following the Priestley-Taylor
    def PLE_plantfun(s, Ac, Fwet, gama):
        PLE = s.expression('(alfa*s*Ac*(1-fwet)/(s+gama))', {
            'alfa': 1.26,
            's': s,
            'Ac': Ac,
            'fwet': Fwet,
            'gama': gama
        })
        return PLE

    # Function to calculate total soil resistance to water vapor transport
    def rtotfun(temp, vpd, bplut, pa):
        ta = temp.expression('pow(((t+273.15)/293.15),1.75)', {'t': temp})
        rcorr = ta.expression('1.0/((101300.0/pa)*ta)', {
            'pa': pa,
            'ta': ta
        })

        cond1 = ee.Image(0)
        cond2 = ee.Image(0)
        cond3 = ee.Image(0)
        cond1x = vpd.lte(bplut.select('VPDopen'))
        cond1 = cond1.where(cond1x, 1)
        cond2x = vpd.gte(bplut.select('VPDclose'))
        cond2 = cond2.where(cond2x, 1)
        cond3x = vpd.gt(bplut.select('VPDopen')).And(vpd.lt(bplut.select('VPDclose')))
        cond3 = cond3.where(cond3x, 1)

        rt1 = bplut.select('RBL_max').multiply(cond1)
        rt2 = bplut.select('RBL_min').multiply(cond2)
        rt3 = bplut.expression('rblmax-((rblmax-rblmin)*(vpdclose-vpd))/(vpdclose-vpdopen)', {
            'rblmax': bplut.select('RBL_max'),
            'rblmin': bplut.select('RBL_min'),
            'vpdclose': bplut.select('VPDclose'),
            'vpd': vpd,
            'vpdopen': bplut.select('VPDopen')
        })

        rt3 = rt3.multiply(cond3)

        rtotc = rt1.add(rt2).add(rt3)
        rtot = rcorr.multiply(rtotc)
        return rtot

    # Function to calculate resistance to radiative heat transfer
    def rrsfun(rho, temp):
        rrs = rho.expression('rho*Cp/(4.0*SBC*T*T*T)', {
            'rho': rho,
            'Cp': Cp,
            'SBC': SBC,
            'T': temp.add(273.15)
        })
        return rrs

    # Function to calculate latent heat flux at saturated soil surface
    def LEwetsoilfun(rtot, rrs, rho, vpd, s, Asoil, Fc, Fwet, gama, daylength):
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
    def PLE_soilfun(rtot, rrs, s, Asoil, rho, Fc, vpd, Fwet, daylength, gama):
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
    def LEsoilfun(rh, vpd, LEwetsoil, PLEsoil):
        para1 = rh.expression('pow((rh*0.01),(vpd/beta))', {
            'rh': rh,
            'vpd': vpd,
            'beta': 250
        })
        LEsoil = PLEsoil.multiply(para1).add(LEwetsoil)
        return LEsoil

    # Aggregate all functions to calculate total ET
    def calc_et(img):
        swrad = img.select('srad')  # W/m2
        albedo = img.select('albedo').multiply(0.001)
        ta = img.select('Tavg')
        ta_day = ta  # degree C
        ta_night = img.select('tmmn').subtract(273.15)
        rhmax = img.select('rmax')
        rhmin = img.select('rmin')  # rh_night
        rh = (rhmax.add(rhmin)).multiply(0.5)  # % rh_day
        vpd_day = img.select('vpd').multiply(1000)  # Pa
        daylength = img.select('dayl')
        nightlength = daylength.expression('24*3600.0-dayl', {'dayl': daylength})
        Fc = img.select('Fc')
        LAI = img.select('LAI')
        svp_night = svpfun(ta_night)  # Pa
        svp_day = svpfun(ta_day)
        vpd_night = svp_night.multiply(ee.Image(1).subtract((rhmin.multiply(0.01))))  # Pa
        pa = Pafun(elev)

        # calculate each parameter values in day and night, respectively
        Rn_day = Rnfun(ta_day, albedo, swrad, daylength)  # J/day/m2
        Rn_night = Rnfun(ta_night, albedo, 0, nightlength)

        Gsoil_day = Gsoilfun(bplut, Rn_day, ta_day, daylength, ta_day, ta_night, Fc)
        Gsoil_night = Gsoilfun(bplut, Rn_night, ta_night, nightlength, ta_day, ta_night, Fc)

        # following the User Guide (2017 version)
        G_day = Gsoil_day.multiply(ee.Image(1).subtract(Fc))
        cond1 = Rn_day.lt(G_day)
        G_day = G_day.where(cond1, Rn_day)

        G_night = Gsoil_night.multiply(ee.Image(1).subtract(Fc))
        cond2 = (Rn_night.subtract(G_night)).lt(Rn_day.multiply(-0.5))
        g1 = Rn_night.add(Rn_day.multiply(0.5))
        G_night = G_night.where(cond2, g1)

        Ac_day = Acfun(Rn_day, Fc)
        Ac_night = Acfun(Rn_night, Fc)
        Asoil_day = Asoilfun(Rn_day, Fc, G_day)
        Asoil_night = Asoilfun(Rn_night, Fc, G_night)

        Fwet_day = Fwetfun(rh)
        Fwet_night = Fwetfun(rhmin)

        rho_day = rhofun(rh, pa, ta_day)
        rho_night = rhofun(rhmin, pa, ta_night)

        lamda_day = lamdafun(ta_day)
        lamda_night = lamdafun(ta_night)

        s_day = sfun(svp_day, ta_day)
        s_night = sfun(svp_night, ta_night)

        LEwet_day = LEwetfun(Ac_day, Fc, rho_day, s_day, vpd_day, daylength, Fwet_day, pa, lamda_day, bplut, LAI,
                             ta_day)
        LEwet_night = LEwetfun(Ac_night, Fc, rho_night, s_night, vpd_night, nightlength, Fwet_night, pa, lamda_night,
                               bplut, LAI, ta_night)

        vpdsc_day = vpd_scalar(bplut, vpd_day)
        vpdsc_night = vpd_scalar(bplut, vpd_night)
        tasc_day = T_scalar(bplut, ta_day)
        tasc_night = T_scalar(bplut, ta_night)
        gama_day = gamafun(pa, lamda_day)
        gama_night = gamafun(pa, lamda_night)

        LEtrans_day = LEtransfun(ta_day, pa, bplut, LAI, Fwet_day, rho_day, lamda_day, gama_day, s_day, Ac_day, Fc,
                                 vpd_day, daylength)
        LEtrans_night = LEtransfun2(ta_night, pa, bplut, LAI, Fwet_night, rho_night, lamda_night, gama_night, s_night,
                                    Ac_night, Fc, vpd_night, nightlength)
        PLEtrans_day = PLE_plantfun(s_day, Ac_day, Fwet_day, gama_day)
        PLEtrans_night = PLE_plantfun(s_night, Ac_night, Fwet_night, gama_night)

        rtot_day = rtotfun(ta_day, vpd_day, bplut, pa)
        rtot_night = rtotfun(ta_night, vpd_night, bplut, pa)
        rrs_day = rrsfun(rho_day, ta_day)
        rrs_night = rrsfun(rho_night, ta_night)

        LEwetsoil_day = LEwetsoilfun(rtot_day, rrs_day, rho_day, vpd_day, s_day, Asoil_day, Fc, Fwet_day, gama_day,
                                     daylength)
        LEwetsoil_night = LEwetsoilfun(rtot_night, rrs_night, rho_night, vpd_night, s_night, Asoil_night, Fc,
                                       Fwet_night, gama_night, nightlength)
        PLEsoil_day = PLE_soilfun(rtot_day, rrs_day, s_day, Asoil_day, rho_day, Fc, vpd_day, Fwet_day, daylength,
                                  gama_day)
        PLEsoil_night = PLE_soilfun(rtot_night, rrs_night, s_night, Asoil_night, rho_night, Fc, vpd_night, Fwet_night,
                                    nightlength, gama_night)
        LEsoil_day = LEsoilfun(rh, vpd_day, LEwetsoil_day, PLEsoil_day)
        LEsoil_night = LEsoilfun(rhmin, vpd_night, LEwetsoil_night, PLEsoil_night)

        # calculate total daily LE, PLE, ET and PET.
        LE_day = LEwet_day.add(LEtrans_day).add(LEsoil_day)
        ET_day = LE_day.divide(lamda_day)
        LE_night = LEwet_night.add(LEtrans_night).add(LEsoil_night)
        ET_night = LE_night.divide(lamda_night)
        LE = LE_day.add(LE_night)  # J/m2
        ET = ET_day.add(ET_night)  # mm/d
        LE_v1 = LE.divide(daylength)  # W/m2

        PLE_day = LEwet_day.add(PLEtrans_day).add(LEwetsoil_day).add(PLEsoil_day)
        PLE_night = LEwet_night.add(PLEtrans_night).add(LEwetsoil_night).add(PLEsoil_night)
        PET_day = PLE_day.divide(lamda_day)
        PET_night = PLE_night.divide(lamda_night)

        doy = ee.Date(img.get('system:time_start'))
        return Rn_day.addBands(Rn_night) \
            .addBands(G_day) \
            .addBands(G_night) \
            .addBands(Ac_day) \
            .addBands(Ac_night) \
            .addBands(Asoil_day) \
            .addBands(Asoil_night) \
            .addBands(Fwet_day) \
            .addBands(Fwet_night) \
            .addBands(rho_day) \
            .addBands(rho_night) \
            .addBands(lamda_day) \
            .addBands(lamda_night) \
            .addBands(s_day) \
            .addBands(s_night) \
            .addBands(LEwet_day) \
            .addBands(LEwet_night) \
            .addBands(gama_day) \
            .addBands(gama_night) \
            .addBands(LEtrans_day) \
            .addBands(LEtrans_night) \
            .addBands(rtot_day) \
            .addBands(rtot_night) \
            .addBands(rrs_day) \
            .addBands(rrs_night) \
            .addBands(LEwetsoil_day) \
            .addBands(LEwetsoil_night) \
            .addBands(LEsoil_day) \
            .addBands(LEsoil_night) \
            .addBands(LE_v1) \
            .addBands(ET) \
            .addBands(daylength) \
            .addBands(LAI) \
            .addBands(Fc) \
            .addBands(ET_day) \
            .addBands(ET_night) \
            .select(
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
             29, 30, 31, 32, 33, 34, 35, 36],
            ['Rn_day', 'Rn_night', 'G_day', 'G_night', 'Ac_day', 'Ac_night', 'Asoil_day', 'Asoil_night', 'Fwet_day',
             'Fwet_night',
             'rho_day', 'rho_night', 'lamda_day', 'lamda_night', 's_day', 's_night', 'LEwet_day', 'LEwet_night',
             'gama_day', 'gama_night',
             'LEtrans_day', 'LEtrans_night', 'rtot_day', 'rtot_night', 'rrs_day', 'rrs_night', 'LEwetsoil_day',
             'LEwetsoil_night',
             'LEsoil_day', 'LEsoil_night', 'LE', 'ET', 'daylength', 'LAI', 'Fc', 'ET_day', 'ET_night']) \
            .copyProperties(img, [time]) \
            .set({'Date': doy})

    return meteo.map(calc_et)
