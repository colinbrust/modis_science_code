import ee
from scripts.MOD16 import MOD16

roi = ee.Geometry.Polygon(
        [[[-112.241, 33.078],
          [-112.241, 31.733],
          [-108.758, 31.733],
          [-108.758, 33.078]]], None, False)
pnt = ee.Geometry.Point([-111.651, 32.631])



year = 2015
meteo = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')

# ----------------------------MOD16 Example-------------------------------------

# Here the Daymet daylength product is used to differentiate between daytime and nighttime ET.
# However, Daymet is only available for North America. In the original MOD16 algorithm, daytime 
# is defined as the period of time when incoming SW radiation is > 10 W/m^2.
dayl = ee.ImageCollection('NASA/ORNL/DAYMET_V3').select('dayl')
elev = ee.Image('USGS/NED')


ET = MOD16(roi, year, meteo, dayl, elev)

print(ET)