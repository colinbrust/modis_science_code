import ee
from msc.utils.GapFill import gap_fill as gf

ee.Initialize()
time = 'system:time_start'
day = (24 * 60 * 60 * 1000)

# Converts 4-day FPAR product to a daily timestep.
# Also mask out problem pixels and temporally gap fill.
def modis_fpar_lai(roi, year):

    start = ee.Date.fromYMD(year,1,1)
    end = ee.Date.fromYMD(year+1,1,1)
    tempwin = 4  # day difference between FPAR images
    tempdif = tempwin*day # Temporal difference in milliseconds

    # Function to filter out cloudy and saturated MODIS FPAR images
    def qMask(img):

        mask1 = img.select('FparLai_QC').bitwiseAnd(32).eq(0)
        mask2 = img.select('FparLai_QC').bitwiseAnd(8).eq(0)
        mask_p = mask1.And(mask2)
        img_masked = img.updateMask(mask_p)

        return img_masked

    # Mask out problem pixels
    fuseEVI = ee.ImageCollection('MODIS/006/MCD15A3H')\
      .filterDate(start.advance(-1, 'month'), end)\
      .map(lambda img:  img.clip(roi))\
      .map(qMask)

    # Fill gaps with temporally nearest valid pixel to create continuous image.
    fuseEVI = gf(fuseEVI)\
    .filterDate(start, end)

    proj = ee.Image(fuseEVI.first()).projection()

    def make_daily_coll(tup):

      return ee.Image.constant(tup).toInt64()\
        .set('system:time_start', ee.Number(tup))\
        .select(['constant'], ['time'])

    # Convert 4-day images to a daily timestep
    ddStart = ee.Number(fuseEVI.aggregate_min('system:time_start'))
    ddEnd =  ee.Number(fuseEVI.aggregate_max('system:time_start'))
    ddList = ee.List.sequence(ddStart, ddEnd, (24*60*60*1000))
    ddCll = ee.ImageCollection.fromImages(ddList.map(lambda tup: make_daily_coll(tup)))

    # function to make collection of days in time period
    def copyValue(img):
        time = img.metadata('system:time_start')

        def mask(val):
          timeOrig = val.metadata('system:time_start')
          masked = timeOrig.eq(time)
          return val.mask(masked)

        eviCll = fuseEVI.map(mask)
        return img.addBands(eviCll.max())

    filledDate = ddCll.map(copyValue)
    # Create a list of original calculated dates
    startDayList = ee.List.sequence(ddStart,ddEnd.subtract(day),tempdif)

    # convert it to an image collection
    def toImage(tup):
      return ee.Image.constant(ee.Number(tup)).set('system:time_start',tup)

    imgCll = ee.ImageCollection.fromImages(startDayList.map(toImage))

    # Function to interpolate 4-day values to daily values
    def ic_interp(ic, band):

        # calculate and interpolate values
        def interpolator(img):
            # get the first day of the subsample
            begin = ee.Number(img.get('system:time_start'))
            # get the end day of the subsample
            end = begin.add(ee.Number(tempdif))
            # convert to an image
            minDD = ee.Image(filledDate.filterDate(begin).first())

            maxDD = ee.Image(filledDate.filterDate(end).first())
            # calculate the coefficent
            angularCoeff = (maxDD.select(band).subtract(minDD.select(band))).divide(tempdif)
            q = ((maxDD.select('time').multiply(minDD.select(band))).subtract(minDD.select('time').multiply(maxDD.select(band)))).divide(tempdif)
            # refilter the out collection
            mnth = filledDate.filterDate(begin,end)

            # interpolate the values
            def ndviInterpolator(img):
              NDVI = (img.select('time').multiply(angularCoeff)).add(q)
              result = img.select(band).unmask(NDVI)

              return result.clip(roi).toShort()

            return mnth.map(ndviInterpolator)

        blockColl = ee.ImageCollection(ic.map(interpolator).flatten().toList(365))

        endDay = blockColl.limit(1, 'system:time_start', False).first()
        end_time = ee.Number(endDay.get('system:time_start'))
        endDay = endDay.set('system:time_start', end_time.add(day))

        return ee.ImageCollection(blockColl.merge(endDay))

    # compute daily interpolations of LAI and Fpar
    LAI = ic_interp(imgCll, 'Lai').map(lambda img: img.rename('LAI').copyProperties(img,['system:time_start']))

    Fc = ic_interp(imgCll, 'Fpar').map(lambda img: img.multiply(0.01).rename('Fc').copyProperties(img, ['system:time_start']))

    return ee.Dictionary({'LAI': LAI,
                          'FPAR': Fc})
