import ee
from msc.calval.EvaluateCalibration import get_optimal_params
from msc.utils.Bplut import make_custom_bplut
from msc.models.GeeMod16 import MOD16
from msc.utils.ExtractGeeData import ExtractGeeData

ee.Initialize()

# Region of interest centered around US-Me2 AmeriFlux site.
roi = ee.Geometry.Point([-121.557, 44.4523]).buffer(100000)

# Bring in MODIS Type 2 landcover classification
ic = ee.ImageCollection('MODIS/006/MCD12Q1') \
        .filterDate('2015-01-01', '2015-12-31') \
        .select('LC_Type2') \
        .first() \
        .clip(roi)

# Calibrated optimal parameters for the model
params = get_optimal_params('data/MOD16/kfold_results',  'parvpd_close > parvpd_open & parrbl_min < parrbl_max')

# Convert calibrated parameters to a dictionary
bp_dict = params.set_index('group').T.to_dict()

# Dictionary that assigns each PFT in bp_dict to the MODIS specified value
mapping_dict = {'ENF': 1, 'EBF': 2, 'DNF': 3, 'DBF': 4, 'MF': 5, 'CSH': 6,
                'OSH': 7, 'WSA': 8, 'SAV': 9, 'GRA': 10, 'CRO': 12}

# Create DNF and SAV from ENF and WSA as there are no tower sites for these PFTs
bp_dict['DNF'] = bp_dict['ENF']
bp_dict['SAV'] = bp_dict['WSA']

# Create the new ee.Image bplut to use as input into the MOD16 model
bplut = make_custom_bplut(ic, bp_dict, mapping_dict)

meteo = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')
dayl = ee.ImageCollection('NASA/ORNL/DAYMET_V3').select('dayl')
elev = ee.Image('USGS/NED')

# Create instance of ExtractGeeData class
extractor = ExtractGeeData(template='data/MOD16/template_example.csv',
                           model=MOD16,
                           out_dir='data/MOD16',
                           meteorology=meteo,
                           daylength=dayl,
                           elev=elev,
                           bplut=bplut)

result = extractor._run_single_location(name='US-AR1', year=2010, lat=36.4267, lon=-99.42)
print(result)
