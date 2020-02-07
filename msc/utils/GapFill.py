import ee
ee.Initialize()


def gap_fill(collection: ee.ImageCollection) -> ee.ImageCollection:
    """
    Function to temporally gap fill images with faulty pixels from QA/QC procedure. Solution taken from
        https://gis.stackexchange.com/questions/293577/sequential-cloud-mask-closest-to-date-in-google-earth-engine.
    :param collection: ee.ImageCollection for temporal gap filling
    :return:ee.ImageCollection with pixels gap filled
    """
    cNew = collection.sort('system:time_start', True)
    collist = cNew.toList(cNew.size())
    first = ee.Image(collist.get(0)).unmask()
    rest = collist.slice(1)

    def wrap(img, ini):
        ini = ee.List(ini)
        img = ee.Image(img)
        last = ee.Image(ini.get(-1))
        mask = img.mask().Not()
        last_masked = last.updateMask(mask)
        last2add = last_masked.unmask()
        img2add = img.unmask()
        added = img2add.add(last2add).set('system:index', ee.String(img.id()))

        props = img.propertyNames()
        condition = props.contains('system:time_start')

        final = ee.Image(ee.Algorithms.If(condition, added.set('system:time_start',
                                                              img.date().millis()), added)).toByte()

        return ini.add(final.copyProperties(img))

    newcol = ee.List(rest.iterate(wrap, ee.List([first])))
    return ee.ImageCollection.fromImages(newcol)
