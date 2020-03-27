import ee
ee.Initialize()


def calc_rh(img, temp_band, band_name):
    rh = img.expression('0.263 * p * q * pow(exp(17.67 * (temp - 273.16)/(temp - 29.65)), -1)', {
        'p': img.select('PS'),
        'q': img.select('QV2M'),
        'temp': img.select(temp_band)
    }).rename(band_name).copyProperties(img, ['system:time_start'])

    rh = ee.Image(rh).clamp(0, 100)
    return img.addBands(rh)


def calc_ravg(img):
    ravg = img.expression('((rmin+rmax)/2)*0.01',{
        'rmin': img.select('rmin'),
        'rmax': img.select('rmax'),
    }).rename('ravg').copyProperties(img, ['system:time_start'])

    ravg = ee.Image(ravg)
    return img.addBands(ravg)


def calc_tavg(img):
    tavg = img.expression('(tmax + tmin)/2', {
        'tmax': img.select('T2MMAX'),
        'tmin': img.select('T2MMIN')
    }).rename('tavg').copyProperties(img, ['system:time_start'])

    tavg = ee.Image(tavg)
    return img.addBands(tavg)


def calc_sat_vp(img, t_band, band_name):
    sat_vp = img.expression('0.6108 * exp((17.27 * t) / (t + 237.3))',{
        't': img.select(t_band)
    }).rename(band_name).copyProperties(img, ['system:time_start'])

    sat_vp = ee.Image(sat_vp)
    return img.addBands(sat_vp)


def calc_vpd(img):
  
    img = calc_sat_vp(img, 'T2MMIN', 'svp_tmin')
    img = calc_sat_vp(img, 'T2MMAX', 'svp_tmax')
    img = calc_sat_vp(img, 'tavg', 'svp_tavg')

    vpa = img.expression('((svp_tmin * rmax/100) + (svp_tmax * rmin/100)) / 2', {
        'svp_tmin': img.select('svp_tmin'),
        'svp_tmax': img.select('svp_tmax'),
        'rmax': img.select('rmax'),
        'rmin': img.select('rmin')
    }).rename('vp_actual').copyProperties(img, ['system:time_start'])

    vpa = ee.Image(vpa)
    img = img.addBands(vpa)

    vpd = img.expression('(svp_tavg - vp_actual)/1000', {
        'svp_tavg': img.select('svp_tavg'),
        'vp_actual': img.select('vp_actual')
    }).rename('vpd').copyProperties(img, ['system:time_start'])

    vpd = ee.Image(vpd)
    return img.addBands(vpd)


def add_all_bands(img):
  
    img = calc_rh(img, 'T2MMAX', 'rmin')
    img = calc_rh(img, 'T2MMIN', 'rmax')
    img = calc_ravg(img)
    img = calc_tavg(img)
    img = calc_vpd(img)

    img = ee.Image(img.copyProperties(img, ['system:time_start']))
    img = img.select(['T2MMAX', 'T2MMIN', 'SWGDN', 'rmin', 'rmax', 'vpd'],
                     ['tmmx', 'tmmn', 'srad', 'rmin', 'rmax', 'vpd'])
    return img


def meteo_proj(img, proj):
    img = img.resample('bicubic').reproject(crs=proj)\
        .copyProperties(img, ['system:time_start', 'system:index'])

    return img


def calc_meteo(ic):

    proj = ic.first().projection()
    out = ic.map(lambda img: meteo_proj(img, proj))

    return out.map(add_all_bands)
