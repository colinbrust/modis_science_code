import ee
ee.Initialize()
import pickle
import pandas as pd
from msc.models.GeeMod17 import MOD17
from msc.models.GeeMod16 import MOD16
from msc.utils.WriteModelToGee import WriteModelToGee
from msc.utils.Smap import get_smap
from msc.calval.EvaluateCalibration import get_optimal_params
from msc.utils.Bplut import make_custom_bplut


def merge_l2_l5_landcover(img):
    lc2 = img.select('LC_Type2')
    lc5 = img.select('LC_Type5').multiply(3)
    out = lc2.where(lc2.eq(12), lc5)

    mask = out.eq(1).Or(out.eq(2)).Or(out.eq(3)).Or(out.eq(4)).Or(out.eq(5)).Or(out.eq(6)).Or(out.eq(7)).Or(out.eq(8)) \
        .Or(out.eq(9)).Or(out.eq(10)).Or(out.eq(21)).Or(out.eq(24))

    return out.updateMask(mask)


# Meteorology required to run MOD17
meteo = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')
dayl = ee.ImageCollection('NASA/ORNL/DAYMET_V3').select('dayl')
elev = ee.Image('USGS/NED')
lc = ee.ImageCollection('MODIS/006/MCD12Q1').map(merge_l2_l5_landcover).filterDate('2015-01-01', '2016-01-01').first()


# Open pickled file of geometries to run the model in
with(open('../data/conus_mosaic.p', 'rb')) as f:
    bounds = pickle.load(f)

# Only use the first 5 geometries to reduce run time
bounds = bounds[0:2]

# Bring in SMAP data
smap = get_smap(start='2003-01-01', end='2018-01-01')

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

bplut = make_custom_bplut(lc, mapping_dict=mapping_dict, bp_dict=bp_dict)

kwargs = {'meteorology': meteo, 'daylength': dayl, 'elev': elev, 'smap_sm': smap, 'bplut': bplut}
# Make instance of WriteModelToGee class
runner = WriteModelToGee(model=MOD16, year=2015, bounds=bounds, unq_id='this_is_a_test',
                         asset_path='users/colinbrust/asdf', **kwargs)

# Run the model, check for progress every 5 minutes
runner.run_model(wait_time=300)
