import ee


# Functions to interpolate data between missing values
# Taken from https://github.com/kongdd/gee_packages/blob/master/Math/pkg_smooth.js and converted to python.
def setweights(ImgCol, bound, ymin):
    if bound is None:
        alpha = 1  # unit: %
        bound = ImgCol.reduce(ee.Reducer.percentile([alpha / 2, 100 - alpha / 2]))

    ymin = ymin or 0

    def mapper(img):
        w = img.multiply(0)
        con_norm = img.expression('b() >= min && b() <= max',
                                  {
                                      min: bound.select(0).max(ymin),  # min should >= ymin
                                      max: bound.select(1)
                                  })
        w = w.where(con_norm, 1)

        return w

    w = ImgCol.map(mapper)
    w = w.toArray()  # .toArray(1)

    return w


def replace_mask(img, newimg, nodata):
    nodata = nodata or 0

    mask = img.mask()

    img = img.unmask(nodata)
    img = img.where(mask.Not(), newimg)

    img = img.updateMask(img.neq(nodata))
    return img


def addTimeBand(img):
    mask = img.mask()
    time = img.metadata('system:time_start').rename("time").mask(mask)
    return img.addBands(time)


def linearInterp(imgcol, frame, nodata):
    frame = frame or 32
    nodata = nodata or 0

    # frame = 32
    time = 'system:time_start'
    imgcol = imgcol.map(addTimeBand)

    # We'll look for all images up to 32 days away from the current image.
    maxDiff = ee.Filter.maxDifference(frame * (1000 * 60 * 60 * 24), time, None, time)

    # Images after, sorted in descending order (so closest is last).
    # f1 = maxDiff.and(ee.Filter.lessThanOrEquals(time, null, time))
    f1 = ee.Filter.And(maxDiff, ee.Filter.lessThanOrEquals(leftField=time, rightField=time))
    c1 = ee.Join.saveAll(matchesKey='after', ordering='time', ascending=False).apply(imgcol, imgcol, f1)

    # Images before, sorted in ascending order (so closest is last).
    # f2 = maxDiff.and(ee.Filter.greaterThanOrEquals(time, null, time))
    f2 = ee.Filter.And(maxDiff, ee.Filter.greaterThanOrEquals(leftField=time, rightField=time))
    c2 = ee.Join.saveAll(matchesKey='before', ordering='time', ascending=True).apply(c1, imgcol, f2)

    def mapper(img):
        img = ee.Image(img)

        before = ee.ImageCollection.fromImages(ee.List(img.get('before'))).mosaic()
        after = ee.ImageCollection.fromImages(ee.List(img.get('after'))).mosaic()

        img = img.set('before', None).set('after', None)
        # constrain after or before no NA values, confirm linear Interp having result
        before = replace_mask(before, after, nodata)
        after = replace_mask(after, before, nodata)

        # Compute the ratio between the image times.
        x1 = before.select('time').double()
        x2 = after.select('time').double()
        now = ee.Image.constant(img.date().millis()).double()
        ratio = now.subtract(x1).divide(x2.subtract(x1))  # this is zero anywhere x1 = x2
        # Compute the interpolated image.
        before = before.select(0)  # remove time band now
        after = after.select(0)
        img = img.select(0)

        interp = after.subtract(before).multiply(ratio).add(before)
        # mask   = img.select([0]).mask()

        qc = img.mask().Not().rename('qc')
        interp = replace_mask(img, interp, nodata)
        # Map.addLayer(interp, {}, 'interp')
        return interp.addBands(qc).copyProperties(img, img.propertyNames())

    interpolated = ee.ImageCollection(c2.map(mapper))

    return interpolated


