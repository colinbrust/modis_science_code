import ee
from scripts import GetModisFpar as fp
from scripts import Bplut as bp

ee.Initialize()


def MOD17(roi, year, meteorology):
    time = 'system:time_start'
    start = ee.Date.fromYMD(year, 1, 1)
    end = ee.Date.fromYMD(year + 1, 1, 1)
    day = (24 * 60 * 60 * 1000)  # Day in milliseconds
    tempwin = 4  # day difference between FPAR images
    tempdif = tempwin * day  # Temporal difference between images in milliseconds

    # Bring in bplut of MOD17 parameters
    bplut = bp.m17_bplut(roi, start, end)

    # Bring in FPAR imagery for the year
    fp_lai = fp.modis_fpar_lai(roi, year)

    fpar = ee.ImageCollection(fp_lai.get('FPAR')).map(lambda img: img.rename('fpar'))
    proj = ee.Image(fpar.first()).projection()

    # Function to downscale inputs to match MODIS projection and resolution
    def match_proj(img):
        img = img.resample('bilinear').reproject(
            crs=proj.crs(),
            scale=500
        ).copyProperties(img, ['system:time_start', 'system:index'])

        return img

    # Reproject, clip and select necessary bands of meteorology
    meteorology = meteorology.filterDate(start, end) \
        .map(lambda img: img.clip(roi)) \
        .map(match_proj) \
        .select(['vpd', 'tmmn', 'tmmx', 'srad'])

    # Function to calculate light use efficiency. Temperature is in K, 
    # pressure is in Pa. 
    def calc_lue(bplut, img):
        Tscalar = img.expression('(tmin-tmmin)/(tmmax-tmmin)', {
            'tmin': img.select('tmmn').subtract(273.15),
            'tmmin': bplut.select('tmmin'),
            'tmmax': bplut.select('tmmax'),
        }).clamp(0, 1)

        VPDscalar = img.expression('(vpdmax-vpd)/(vpdmax-vpdmin)', {
            'vpd': img.select('vpd').multiply(1000),
            'vpdmax': bplut.select('vpdmax'),
            'vpdmin': bplut.select('vpdmin'),
        }).clamp(0, 1)

        lue = bplut.expression('luemax*Tscalar*VPDscalar', {
            'luemax': bplut.select('luemax'),
            'Tscalar': Tscalar,
            'VPDscalar': VPDscalar,
        })

        return lue.addBands(img.select('srad')).rename(['lue', 'srad'])

    # Join meteorology and lue data
    def data_join(left, right):
        filter = ee.Filter.maxDifference(
            difference=tempdif,
            leftField=time,
            rightField=time
        )

        join = ee.Join.saveBest(
            matchKey='bestImage',
            measureKey='delta_t'
        ).apply(
            primary=left,
            secondary=right,
            condition=filter
        )

        return ee.ImageCollection(join).map(lambda img: ee.Image(img).addBands(img.get('bestImage')))

    # Function to calculate GPP
    def calc_gpp(gpp_components):
        # Convert incoming SW rad from Watts to MJ to MJ
        par = gpp_components.expression('(srad*0.45*60*60*24*0.000001)', {
            'srad': gpp_components.select('srad'),
        })

        # Calcuate GPP as product of FPAR, PAR, and scaled light use efficiency.
        gpp = gpp_components.expression('(fpar*par*lue)', {
            'fpar': gpp_components.select('fpar'),
            'par': par,
            'lue': gpp_components.select('lue')
        })

        return gpp.rename('GPP')

    # Calculate light use efficiency value given meteorology
    lue_srad = meteorology.map(lambda img: calc_lue(bplut, img).copyProperties(img, [time]))

    # Join data and calculate GPP
    gpp_components = data_join(fpar, lue_srad)
    gpp = gpp_components.map(lambda img: calc_gpp(img).copyProperties(img, [time]))

    gpp_out = data_join(gpp_components, gpp)
    gpp_out = data_join(gpp_out, meteorology)

    return gpp_out
