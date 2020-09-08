import ee
from typing import Dict


def make_custom_bplut(img: ee.Image,
                      bp_dict: Dict[str, Dict[str, float]],
                      mapping_dict: Dict[str, int]) -> ee.Image:
    params = list(list(bp_dict.values())[0].keys())
    pfts = list(bp_dict.keys())

    dict_out = {}

    for param in params:
        fill_dict = {}
        for pft in pfts:
            fill_dict[mapping_dict[pft]] = bp_dict[pft][param]
        dict_out[param] = fill_dict

    spatial_bp = dict(zip(params, [ee.Image(0)] * len(params)))

    for pft in mapping_dict.values():
        for param in spatial_bp.keys():
            spatial_bp[param] = spatial_bp[param].where(img.eq(pft), dict_out[param][pft])

    bp_out = img.rename('landcover')

    for param, img2 in spatial_bp.items():
        bp_out = bp_out.addBands(img2.rename(param))

    return ee.Image(bp_out.updateMask(img.neq(0))).copyProperties(img, ['system:time_start'])


def make_pft_lc(roi, start, end):
    pft = ee.ImageCollection('MODIS/006/MCD12Q1') \
        .filterDate(start, end.advance(-1, 'day')) \
        .select('LC_Type2') \
        .first() \
        .clip(roi)

    # 1 = ENF, 2 = EBF, 3 = DNF, 4 = DBF, 5 = MF, 6 = CSH,
    # 7 = OSH, 8 = WSA, 9 = SAV, 10 = GRA, 11 = CRO

    # Distinguish different MODIS plant functional types
    pft_lc = ee.Image(0)
    pft_lc = pft_lc.where(pft.eq(1), 1)
    pft_lc = pft_lc.where(pft.eq(2), 2)
    pft_lc = pft_lc.where(pft.eq(3), 3)
    pft_lc = pft_lc.where(pft.eq(4), 4)
    pft_lc = pft_lc.where(pft.eq(5), 5)
    pft_lc = pft_lc.where(pft.eq(6), 6)
    pft_lc = pft_lc.where(pft.eq(7), 7)
    pft_lc = pft_lc.where(pft.eq(8), 8)
    pft_lc = pft_lc.where(pft.eq(9), 9)
    pft_lc = pft_lc.where(pft.eq(10), 10)
    pft_lc = pft_lc.where(pft.eq(12), 11)

    return pft_lc


def m16_bplut(roi, start, end):
    pft_lc = make_pft_lc(roi, start, end)

    # Update mask so that we only use images with a valid land cover type
    lc = pft_lc.updateMask(pft_lc.neq(0))

    Tminclose = {1: -8.0,
                 2: -8.0,
                 3: -8.0,
                 4: -6.0,
                 5: -7.0,
                 6: -8.0,
                 7: -8.0,
                 8: -8.0,
                 9: -8.0,
                 10: -8.0,
                 11: -8.0}

    Tminopen = {1: 8.31,
                2: 9.09,
                3: 10.44,
                4: 9.94,
                5: 9.50,
                6: 8.61,
                7: 8.80,
                8: 11.39,
                9: 11.39,
                10: 12.02,
                11: 12.02}

    VPDclose = {1: 3000.0,
                2: 4000.0,
                3: 3500.0,
                4: 2900.0,
                5: 2900.0,
                6: 4300.0,
                7: 4400.0,
                8: 3500.0,
                9: 3600.0,
                10: 4200.0,
                11: 4500.0}

    VPDopen = {1: 650.0,
               2: 1000.0,
               3: 650.0,
               4: 650.0,
               5: 650.0,
               6: 650.0,
               7: 650.0,
               8: 650.0,
               9: 650.0,
               10: 650.0,
               11: 650.0}

    Gl_sh = {1: 0.04,
             2: 0.01,
             3: 0.04,
             4: 0.01,
             5: 0.04,
             6: 0.04,
             7: 0.04,
             8: 0.08,
             9: 0.08,
             10: 0.02,
             11: 0.02}

    Gl_e_wv = {1: 0.04,
               2: 0.01,
               3: 0.04,
               4: 0.01,
               5: 0.04,
               6: 0.04,
               7: 0.04,
               8: 0.08,
               9: 0.08,
               10: 0.02,
               11: 0.02}

    CL = {1: 0.0032,
          2: 0.0025,
          3: 0.0032,
          4: 0.0028,
          5: 0.0025,
          6: 0.0065,
          7: 0.0065,
          8: 0.0065,
          9: 0.0065,
          10: 0.0070,
          11: 0.0070}

    RBL_min = {1: 65.0,
               2: 70.0,
               3: 65.0,
               4: 65.0,
               5: 65.0,
               6: 20.0,
               7: 20.0,
               8: 25.0,
               9: 25.0,
               10: 20.0,
               11: 20.0}

    RBL_max = {1: 95.0,
               2: 100.0,
               3: 95.0,
               4: 100.0,
               5: 95.0,
               6: 55.0,
               7: 55.0,
               8: 45.0,
               9: 45.0,
               10: 50.0,
               11: 50.0}

    SM_open = {1: 0.192,
               2: 0.447,
               3: 0.192,
               4: 0.177,
               5: 0.338,
               6: 0.164,
               7: 0.142,
               8: 0.481,
               9: 0.481,
               10: 0.245,
               11: 0.250}

    SM_close = {1: 0.083,
                2: 0.073,
                3: 0.083,
                4: 0.020,
                5: 0.051,
                6: 0.016,
                7: 0.108,
                8: 0.146,
                9: 0.146,
                10: 0.019,
                11: 0.050}

    Beta = {1: 250.0,
            2: 250.0,
            3: 250.0,
            4: 250.0,
            5: 250.0,
            6: 250.0,
            7: 250.0,
            8: 250.0,
            9: 250.0,
            10: 250.0,
            11: 250.0}

    tminclose = ee.Image(0)
    tminopen = ee.Image(0)
    vpdclose = ee.Image(0)
    vpdopen = ee.Image(0)
    smclose = ee.Image(0)
    smopen = ee.Image(0)
    gl_sh = ee.Image(0)
    gl_e_wv = ee.Image(0)
    cl = ee.Image(0)
    rblmin = ee.Image(0)
    rblmax = ee.Image(0)
    beta = ee.Image(0)

    for i in range(1, 12):
        tminclose = tminclose.where(lc.eq(i), Tminclose[i])
        tminopen = tminopen.where(lc.eq(i), Tminopen[i])
        vpdclose = vpdclose.where(lc.eq(i), VPDclose[i])
        vpdopen = vpdopen.where(lc.eq(i), VPDopen[i])
        smclose = smclose.where(lc.eq(i), SM_close[i])
        smopen = smopen.where(lc.eq(i), SM_open[i])
        gl_sh = gl_sh.where(lc.eq(i), Gl_sh[i])
        gl_e_wv = gl_e_wv.where(lc.eq(i), Gl_e_wv[i])
        cl = cl.where(lc.eq(i), CL[i])
        rblmin = rblmin.where(lc.eq(i), RBL_min[i])
        rblmax = rblmax.where(lc.eq(i), RBL_max[i])
        beta = beta.where(lc.eq(i), Beta[i])

    bplut = lc.addBands(tminclose) \
        .addBands(tminopen) \
        .addBands(vpdclose) \
        .addBands(vpdopen) \
        .addBands(smopen) \
        .addBands(smclose) \
        .addBands(gl_sh) \
        .addBands(gl_e_wv) \
        .addBands(cl) \
        .addBands(rblmin) \
        .addBands(rblmax) \
        .addBands(beta) \
        .select([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                ['landcover', 'tmin_close', 'tmin_open', 'vpd_close', 'vpd_open', 'sm_open', 'sm_close',
                 'gl_sh', 'gl_e_wv', 'Cl', 'rbl_min', 'rbl_max', 'beta'])

    return bplut.updateMask(pft_lc.neq(0))


def zhang_m16_bplut(roi, start, end):
    pft_lc = make_pft_lc(roi, start, end)

    # Update mask so that we only use images with a valid land cover type
    lc = pft_lc.updateMask(pft_lc.neq(0))

    Tminclose = {1: -8.0,
                 2: -8.0,
                 3: -8.0,
                 4: -6.0,
                 5: -7.0,
                 6: -8.0,
                 7: -8.0,
                 8: -8.0,
                 9: -8.0,
                 10: -8.0,
                 11: -8.0}

    Tminopen = {1: 8.31,
                2: 9.09,
                3: 10.44,
                4: 9.94,
                5: 9.50,
                6: 8.61,
                7: 8.80,
                8: 11.39,
                9: 11.39,
                10: 12.02,
                11: 12.02}

    VPDclose = {1: 5100.0,
                2: 4800.0,
                3: 5400.0,
                4: 6000.0,
                5: 5200.0,
                6: 4000.0,
                7: 5700.0,
                8: 5600.0,
                9: 5600.0,
                10: 5800.0,
                11: 5200.0}

    VPDopen = {1: 800.0,
               2: 700.0,
               3: 700.0,
               4: 900.0,
               5: 900.0,
               6: 800.0,
               7: 1300.0,
               8: 1100.0,
               9: 1100.0,
               10: 1200.0,
               11: 1100.0}

    Gl_sh = {1: 0.03,
             2: 0.05,
             3: 0.01,
             4: 0.01,
             5: 0.05,
             6: 0.04,
             7: 0.02,
             8: 0.03,
             9: 0.03,
             10: 0.01,
             11: 0.01}

    Gl_e_wv = {1: 0.007,
               2: 0.006,
               3: 0.013,
               4: 0.006,
               5: 0.007,
               6: 0.006,
               7: 0.007,
               8: 0.009,
               9: 0.009,
               10: 0.006,
               11: 0.007}

    CL = {1: 0.0014,
          2: 0.0014,
          3: 0.0018,
          4: 0.0013,
          5: 0.0014,
          6: 0.0017,
          7: 0.0019,
          8: 0.0075,
          9: 0.0074,
          10: 0.0066,
          11: 0.0074}

    RBL_min = {1: 50.0,
               2: 55.0,
               3: 40.0,
               4: 20.0,
               5: 55.0,
               6: 65.0,
               7: 45.0,
               8: 60.0,
               9: 60.0,
               10: 70.0,
               11: 50.0}

    RBL_max = {1: 110.0,
               2: 75.0,
               3: 140.0,
               4: 100.0,
               5: 120.0,
               6: 90.0,
               7: 140.0,
               8: 115.0,
               9: 115.0,
               10: 135.0,
               11: 130.0}

    Beta = {1: 900.0,
            2: 800.0,
            3: 300.0,
            4: 900.0,
            5: 800.0,
            6: 400.0,
            7: 500.0,
            8: 500.0,
            9: 600.0,
            10: 800.0,
            11: 900.0}

    SM_open = {1: 0.192,
               2: 0.447,
               3: 0.192,
               4: 0.177,
               5: 0.338,
               6: 0.164,
               7: 0.142,
               8: 0.481,
               9: 0.481,
               10: 0.245,
               11: 0.250}

    SM_close = {1: 0.083,
                2: 0.073,
                3: 0.083,
                4: 0.020,
                5: 0.051,
                6: 0.016,
                7: 0.108,
                8: 0.146,
                9: 0.146,
                10: 0.019,
                11: 0.050}

    tminclose = ee.Image(0)
    tminopen = ee.Image(0)
    vpdclose = ee.Image(0)
    vpdopen = ee.Image(0)
    smclose = ee.Image(0)
    smopen = ee.Image(0)
    gl_sh = ee.Image(0)
    gl_e_wv = ee.Image(0)
    cl = ee.Image(0)
    rblmin = ee.Image(0)
    rblmax = ee.Image(0)
    beta = ee.Image(0)

    for i in range(1, 12):
        tminclose = tminclose.where(lc.eq(i), Tminclose[i])
        tminopen = tminopen.where(lc.eq(i), Tminopen[i])
        vpdclose = vpdclose.where(lc.eq(i), VPDclose[i])
        vpdopen = vpdopen.where(lc.eq(i), VPDopen[i])
        smclose = smclose.where(lc.eq(i), SM_close[i])
        smopen = smopen.where(lc.eq(i), SM_open[i])
        gl_sh = gl_sh.where(lc.eq(i), Gl_sh[i])
        gl_e_wv = gl_e_wv.where(lc.eq(i), Gl_e_wv[i])
        cl = cl.where(lc.eq(i), CL[i])
        rblmin = rblmin.where(lc.eq(i), RBL_min[i])
        rblmax = rblmax.where(lc.eq(i), RBL_max[i])
        beta = beta.where(lc.eq(i), Beta[i])

    bplut = lc.addBands(tminclose) \
        .addBands(tminopen) \
        .addBands(vpdclose) \
        .addBands(vpdopen) \
        .addBands(smopen) \
        .addBands(smclose) \
        .addBands(gl_sh) \
        .addBands(gl_e_wv) \
        .addBands(cl) \
        .addBands(rblmin) \
        .addBands(rblmax) \
        .addBands(beta) \
        .select([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                ['landcover', 'tmin_close', 'tmin_open', 'vpd_close', 'vpd_open', 'sm_open', 'sm_close',
                 'gl_sh', 'gl_e_wv', 'Cl', 'rbl_min', 'rbl_max', 'beta'])

    return bplut.updateMask(pft_lc.neq(0))


def m17_bplut(roi, start, end):
    pft_lc = make_pft_lc(roi, start, end)

    # Update mask so that we only use images with a valid land cover type
    lc = pft_lc.updateMask(pft_lc.neq(0))

    LUEmax = {1: 0.000962,
              2: 0.001268,
              3: 0.001086,
              4: 0.001165,
              5: 0.001051,
              6: 0.001281,
              7: 0.000841,
              8: 0.001239,
              9: 0.001206,
              10: 0.000860,
              11: 0.001044}

    Tmin_min = {1: -8.0,
                2: -8.0,
                3: -8.0,
                4: -6.0,
                5: -7.0,
                6: -8.0,
                7: -8.0,
                8: -8.0,
                9: -8.0,
                10: -8.0,
                11: -8.0}

    Tmin_max = {1: 8.31,
                2: 9.09,
                3: 10.44,
                4: 9.94,
                5: 9.50,
                6: 8.61,
                7: 8.80,
                8: 11.39,
                9: 11.39,
                10: 12.02,
                11: 12.02}

    VPD_min = {1: 650.0,
               2: 800.0,
               3: 650.0,
               4: 650.0,
               5: 650.0,
               6: 650.0,
               7: 650.0,
               8: 650.0,
               9: 650.0,
               10: 650.0,
               11: 650.0}

    VPD_max = {1: 4600.0,
               2: 3100.0,
               3: 2300.0,
               4: 1650.0,
               5: 2400.0,
               6: 4700.0,
               7: 4800.0,
               8: 3200.0,
               9: 3100.0,
               10: 5300.0,
               11: 4300.0}

    luemax = ee.Image(0)
    tminmin = ee.Image(0)
    tminmax = ee.Image(0)
    vpdmin = ee.Image(0)
    vpdmax = ee.Image(0)

    for i in range(1, 12):
        luemax = luemax.where(lc.eq(i), LUEmax[i])
        tminmin = tminmin.where(lc.eq(i), Tmin_min[i])
        tminmax = tminmax.where(lc.eq(i), Tmin_max[i])
        vpdmin = vpdmin.where(lc.eq(i), VPD_min[i])
        vpdmax = vpdmax.where(lc.eq(i), VPD_max[i])

    bplut = lc.addBands(luemax) \
        .addBands(tminmin) \
        .addBands(tminmax) \
        .addBands(vpdmin) \
        .addBands(vpdmax) \
        .select([0, 1, 2, 3, 4, 5], ['landcover', 'luemax', 'tmmin', 'tmmax', 'vpdmax', 'vpdmin'])

    return bplut.updateMask(pft_lc.neq(0))
