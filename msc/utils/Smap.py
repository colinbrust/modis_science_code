import ee

ee.Initialize()

sm2015 = ee.ImageCollection('users/colinbrust/l4rz_2015')
sm2016 = ee.ImageCollection('users/colinbrust/l4rz_2016')
sm2017 = ee.ImageCollection('users/colinbrust/l4rz_2017')
sm2018 = ee.ImageCollection('users/colinbrust/l4rz_2018')

nr2000 = ee.ImageCollection('users/colinbrust/nature-run-2000')
nr2001 = ee.ImageCollection('users/colinbrust/nature-run-2001')
nr2002 = ee.ImageCollection('users/colinbrust/nature-run-2002')
nr2003 = ee.ImageCollection('users/colinbrust/nature-run-2003')
nr2004 = ee.ImageCollection('users/colinbrust/nature-run-2004')
nr2005 = ee.ImageCollection('users/colinbrust/nature-run-2005')
nr2006 = ee.ImageCollection('users/colinbrust/nature-run-2006')
nr2007 = ee.ImageCollection('users/colinbrust/nature-run-2007')
nr2008 = ee.ImageCollection('users/colinbrust/nature-run-2008')
nr2009 = ee.ImageCollection('users/colinbrust/nature-run-2009')
nr2010 = ee.ImageCollection('users/colinbrust/nature-run-2010')
nr2011 = ee.ImageCollection('users/colinbrust/nature-run-2011')
nr2012 = ee.ImageCollection('users/colinbrust/nature-run-2012')
nr2013 = ee.ImageCollection('users/colinbrust/nature-run-2013')
nr2014 = ee.ImageCollection('users/colinbrust/nature-run-2014')
nr2015 = ee.ImageCollection('users/colinbrust/nature-run-2015')
nr2016 = ee.ImageCollection('users/colinbrust/nature-run-2016')
nr2017 = ee.ImageCollection('users/colinbrust/nature-run-2017')


def nrMerged():
    return nr2000.merge(nr2001).merge(nr2002).merge(nr2003).merge(nr2004) \
        .merge(nr2005).merge(nr2006).merge(nr2006).merge(nr2007).merge(nr2008) \
        .merge(nr2009).merge(nr2010).merge(nr2011).merge(nr2012).merge(nr2013) \
        .merge(nr2014).merge(nr2015).merge(nr2016).merge(nr2017) \
        .select(['b1', 'b2'],
                ['surfMean', 'rzMean'])


def merged():
    nr_in = nr2015.filterDate('2015-01-01', '2015-04-01') \
        .select(['b1', 'b2'],
                ['surfMean', 'rzMean'])

    return nr_in.merge(sm2015).merge(sm2016).merge(sm2017).merge(sm2018)


def fSM_calc(img, surf_min, surf_max):
    fSM1 = img.expression('(sm - smMin)/(smMax - smMin)', {
        'sm': img.select('surfMean'),
        'smMin': surf_min,
        'smMax': surf_max})

    return fSM1.rename('fSM').copyProperties(img, ['system:index', 'system:time_start'])


def get_smap(start, end, roi):

    smap = nrMerged().map(lambda img: img.clip(roi))

    surf_min = smap.select('surfMean').min()
    surf_max = smap.select('surfMean').max()

    smRz = smap.filterDate(start, end).select('rzMean')
    fSM = smap.filterDate(start, end).map(lambda x: fSM_calc(x, surf_min, surf_max))

    return [smRz, fSM]
