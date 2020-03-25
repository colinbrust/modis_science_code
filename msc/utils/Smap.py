import ee

l4 = ee.ImageCollection('users/colinbrust/L4SM')
nr = ee.ImageCollection('users/colinbrust/NatureRun').filterDate('2009-01-01', '2015-03-31')

SMAP = nr.merge(l4)


def fSM_calc(img, surf_min, surf_max):
    fSM1 = img.expression('(sm - smMin)/(smMax - smMin)', {
        'sm': img.select('sm-surface'),
        'smMin': surf_min,
        'smMax': surf_max}).rename('fSM')

    return img.addBands(fSM1).copyProperties(img, ['system:index', 'system:time_start'])


def get_smap(start, end):

    surf_min = SMAP.select('sm-surface').min()
    surf_max = SMAP.select('sm-surface').max()

    fSM = SMAP.filterDate(start, end).map(lambda x: fSM_calc(x, surf_min, surf_max))

    return fSM
