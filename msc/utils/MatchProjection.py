import argparse
import os
from osgeo import gdal
from pathlib import Path


"""
Taken from: https://github.com/jgomezdans/eoldas_ng_observations/blob/master/eoldas_ng_observations/eoldas_observation_helpers.py#L29
"""


def reproject_image_to_master(master, slave, out_name):
    slave_ds = gdal.Open(slave)
    if slave_ds is None:
        raise IOError("GDAL could not open slave file {}".format(slave))
    slave_proj = slave_ds.GetProjection()
    data_type = slave_ds.GetRasterBand(1).DataType
    n_bands = slave_ds.RasterCount

    master_ds = gdal.Open(master)
    if master_ds is None:
        raise IOError("GDAL could not open master file {}".format(master))
    master_proj = master_ds.GetProjection()
    master_geotrans = master_ds.GetGeoTransform()
    w = master_ds.RasterXSize
    h = master_ds.RasterYSize

    print(out_name)
    dst_ds = gdal.GetDriverByName('GTiff').Create(out_name, w, h, n_bands, data_type)
    dst_ds.SetGeoTransform(master_geotrans)
    dst_ds.SetProjection(master_proj)
    gdal.ReprojectImage(slave_ds, dst_ds, slave_proj, master_proj, gdal.GRA_CubicSpline)
    dst_ds = None  # Flush to disk
    return out_name


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Reproject/Resize a directory of GTiffs to a template proj/size')
    parser.add_argument('template', metavar='t', type=str, nargs=1,
                        help='Template GTiff file')
    parser.add_argument('old_dir', metavar='od', type=str, nargs=1,
                        help='Directory containing images for reprojection')
    parser.add_argument('new_dir', metavar='nd', type=str, nargs=1,
                        help='Directory where new images will be saved')

    args = parser.parse_args()

    for f in Path(args.old_dir[0]).glob('*.tif'):
        new_name = os.path.join(args.new_dir[0], os.path.basename(f))
        reproject_image_to_master(args.template[0], str(f), new_name)


