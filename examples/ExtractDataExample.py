import ee
from scripts.GeeMod16 import MOD16
from scripts.GeeMod17 import MOD17
from scripts.ExtractGeeData import ExtractGeeData

ee.Initialize()

# ---------------------------- MOD16 Example ----------------------------

# # Get ImageCollections necessary to run MOD16 algorithm in GEE
# meteo = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')
# dayl = ee.ImageCollection('NASA/ORNL/DAYMET_V3').select('dayl')
# elev = ee.Image('USGS/NED')
#
# # Create instance of ExtractGeeData class
# extractor = ExtractGeeData(template='data/MOD16/template_example.csv',
#                            model=MOD16,
#                            out_dir='data/MOD16',
#                            meteorology=meteo,
#                            daylength=dayl,
#                            elev=elev)
#
# # Extract data for sites specified in template file
# extractor.run_model(restart=False)

# ---------------------------- MOD17 Example ----------------------------

# Get ImageCollections necessary to run MOD16 algorithm in GEE
meteo = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')


# Create instance of ExtractGeeData class
extractor = ExtractGeeData(template='data/MOD17/template_example.csv',
                           model=MOD17,
                           out_dir='data/MOD17',
                           meteorology=meteo)

# Extract data for sites specified in template file
extractor.run_model(restart=False)
