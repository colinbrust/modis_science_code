import ee

ee.Initialize()
import pandas as pd
from msc.models.GeeMod16 import MOD16
from msc.utils.Smap import get_smap_l4
from msc.calval.EvaluateCalibration import get_optimal_params
from msc.utils.Bplut import make_custom_bplut

start = '2015-01-01'
end = '2016-01-01'


def merge_l2_l5_landcover(img):
    lc2 = img.select('LC_Type2')
    lc5 = img.select('LC_Type5').multiply(3)
    out = lc2.where(lc2.eq(12), lc5)

    mask = out.eq(1).Or(out.eq(2)).Or(out.eq(3)).Or(out.eq(4)).Or(out.eq(5)).Or(out.eq(6)).Or(out.eq(7)).Or(out.eq(8)) \
        .Or(out.eq(9)).Or(out.eq(10)).Or(out.eq(21)).Or(out.eq(24))

    return out.updateMask(mask)


dayl = ee.ImageCollection('NASA/ORNL/DAYMET_V3').select('dayl').filterDate(start, end)
elev = ee.Image('USGS/NED')

gm = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET') \
    .select(['rmax', 'rmin', 'srad', 'tmmx', 'tmmn', 'vpd']) \
    .filterDate(start, end)

fpar = ee.ImageCollection('MODIS/006/MCD15A3H') \
    .select(['Fpar', 'Lai']) \
    .filterDate(start, end)

lc = ee.ImageCollection('MODIS/006/MCD12Q1') \
    .map(merge_l2_l5_landcover) \
    .filterDate(start, end)

smap = get_smap_l4(start=start, end=end)
fc = fpar.select('Fpar')
lai = fpar.select('Lai')

params = get_optimal_params('/mnt/e/Data/GEE_ET/kfold_results/sm_new_format',
                            'parvpd_close > parvpd_open & parrbl_min < parrbl_max & parsm_open > parsm_close')

# Add in default tmin_open and tmin_close params that were not used for calibration
default = pd.read_csv('/mnt/e/PycharmProjects/modis_science_code/data/MOD16/default_bplut.csv')
cols_to_use = [x for x in default.columns.values if x not in params.columns.values] + ['group']
default = default[cols_to_use]

params = pd.merge(params, default, left_on='group', right_on='group')

# Convert calibrated parameters to a dictionary
bp_dict = params.set_index('group').T.to_dict()

bp_dict['DNF'] = bp_dict['ENF']
bp_dict['SAV'] = bp_dict['WSA']

# Dictionary that assigns each PFT in bp_dict to the MODIS specified value
mapping_dict = {'ENF': 1, 'EBF': 2, 'DNF': 3, 'DBF': 4, 'MF': 5, 'CSH': 6,
                'OSH': 7, 'WSA': 8, 'SAV': 9, 'GRA': 10, 'C3': 21, 'C4': 24}

bplut = ee.ImageCollection(lc.map(lambda img: make_custom_bplut(img, mapping_dict=mapping_dict, bp_dict=bp_dict)))

datasets = {'meteorology': gm,
            'bplut': bplut.first(),
            'smapsm': smap,
            'Fc': fc,
            'LAI': lai,
            'daylength': dayl,
            'elev': elev}

roi = ee.Geometry.Polygon(
        [[[-111.01025273380935, 32.00617101764716],
          [-111.01025273380935, 31.485319895486562],
          [-109.48864628849685, 31.485319895486562],
          [-109.48864628849685, 32.00617101764716]]], None, False)

out = MOD16(roi, 2015, **datasets)

for i in range(1, 366):
    day_out = out.filter(ee.Filter.calendarRange(i, i, 'day_of_year')).first()
    task = ee.batch.Export.image.toDrive(
        image=day_out,
        maxPixels=1e13,
        description='SMAP_MOD16_ET_'+str(i),
        fileNamePrefix='SMAP_MOD16_ET_'+str(i),
        scale=500,
        folder='GEE_Exports/SMAPVEX',
        region=[[[-111.01025273380935, 32.00617101764716],
                 [-111.01025273380935, 31.485319895486562],
                 [-109.48864628849685, 31.485319895486562],
                 [-109.48864628849685, 32.00617101764716]]]
    )
    print('SMAP_MOD16_ET_'+str(i))
    task.start()