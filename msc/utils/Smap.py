import ee

l4 = ee.ImageCollection('users/colinbrust/L4SM')
nr = ee.ImageCollection('users/colinbrust/NatureRun').filterDate('2003-01-01', '2015-03-31')

SMAP = nr.merge(l4)

MIN = SMAP.select('sm-surface').min()
MAX = SMAP.select('sm-surface').max()


def fSM_calc(img, surf_min, surf_max):
    fSM1 = img.expression('(sm - smMin)/(smMax - smMin)', {
        'sm': img.select('sm-surface'),
        'smMin': surf_min,
        'smMax': surf_max}).rename('fSM')
    return img.addBands(fSM1).copyProperties(img, ['system:index', 'system:time_start'])


def get_smap(start, end):

    fSM = SMAP.filterDate(start, end).map(lambda x: fSM_calc(x, MIN, MAX))
    return fSM


def get_smap_premade(coll):

    fSM = coll.map(lambda x: fSM_calc(x, MIN, MAX))
    return fSM
