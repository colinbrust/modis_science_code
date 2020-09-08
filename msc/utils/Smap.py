import ee

l4 = ee.ImageCollection('users/colinbrust/L4SM')
nr = ee.ImageCollection('users/colinbrust/NatureRun').filterDate('2003-01-01', '2015-03-31')

SMAP = nr.merge(l4)

SURF_MIN = SMAP.select('sm-surface').min()
SURF_MAX = SMAP.select('sm-surface').max()

ROOT_MIN = SMAP.select('sm-rootzone').min()
ROOT_MAX = SMAP.select('sm-rootzone').max()


def fSM_calc(img, surf_min, surf_max):
    fSM1 = img.expression('(sm - smMin)/(smMax - smMin)', {
        'sm': img.select('sm-surface'),
        'smMin': surf_min,
        'smMax': surf_max}).rename('fSM')
    return img.addBands(fSM1).copyProperties(img, ['system:index', 'system:time_start'])


def fSM_rz_calc(img, rz_min, rz_max):

    fSM1 = img.expression('(sm - smMin)/(smMax - smMin)', {
        'sm': img.select('sm-rootzone'),
        'smMin': rz_min,
        'smMax': rz_max}).rename('fSM-rz')

    return img.addBands(fSM1).copyProperties(img, ['system:index', 'system:time_start'])


def get_smap(start, end):

    fSM = SMAP.filterDate(start, end).map(lambda x: fSM_calc(x, SURF_MIN, SURF_MAX))
    fSM = fSM.map(lambda x: fSM_rz_calc(x, ROOT_MIN, ROOT_MAX))
    return fSM


def get_smap_premade(coll):

    fSM = coll.map(lambda x: fSM_calc(x, SURF_MIN, SURF_MAX))
    fSM = fSM.map(lambda x: fSM_rz_calc(x, ROOT_MIN, ROOT_MAX))

    return fSM
