import ee
from scripts.ee_MOD16 import MOD16
from scripts.ExtractGeeData import ExtractGeeData

year = 2015
meteo = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')

# ----------------------------MOD16 Example-------------------------------------

# Here the Daymet daylength product is used to differentiate between daytime and nighttime ET.
# However, Daymet is only available for North America. In the original MOD16 algorithm, daytime 
# is defined as the period of time when incoming SW radiation is > 10 W/m^2.
dayl = ee.ImageCollection('NASA/ORNL/DAYMET_V3').select('dayl')
elev = ee.Image('USGS/NED')

extractor = ExtractGeeData(template='../data/template_example.csv',
                           model=MOD16,
                           out_dir='../data',
                           meteorology=meteo,
                           daylength=dayl,
                           elev=elev)

extractor.run_model(restart=True)
