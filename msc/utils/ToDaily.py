import ee
from msc.utils.InterpMissing import linearInterp


def to_daily(ic, band, new_band_name):
    new_band_name = new_band_name or band

    ddStart = ee.Number(ic.aggregate_min('system:time_start'))
    ddEnd = ee.Number(ic.aggregate_max('system:time_start'))
    ddList = ee.List.sequence(ddStart, ddEnd, (24 * 60 * 60 * 1000))

    def mapper(tup):
        return ee.Image.constant(0).set('system:time_start', ee.Number(tup)).select(['constant'], ['time'])

    ddCll = ee.ImageCollection.fromImages(ddList.map(mapper))
    time = 'system:time_start'

    def dataJoin(left, right):
        filt = ee.Filter.maxDifference(
            difference=4 * (24 * 60 * 60 * 1000),
            leftField=time,
            rightField=time
        )

        joined = ee.Join.saveBest(
            matchKey='match',
            measureKey='dt',
            outer=True
        )

        return ee.ImageCollection(joined.apply(left, right, filt)) \
            .map(lambda img: img.addBands(img.get('match')).select(band).rename(new_band_name))

    daily = dataJoin(ddCll, ic)

    return linearInterp(daily, 16 * 3, 0).map(lambda img: img.select(new_band_name))
