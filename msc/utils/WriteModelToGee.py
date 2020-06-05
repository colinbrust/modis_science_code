import ee
import os
import re
import shutil
from typing import List, Callable, Dict
import subprocess as sp
import time
from msc.utils.CLIUtils import task_running


class WriteModelToGee(object):
    """
    Class that runs a model across a list of tiles, then mosaics the tiles and aggregates the model results.
    """
    def __init__(self, model: Callable, year: int, bounds: List, unq_id: str, asset_path: str,
                 scale: int, **kwargs: Dict) -> None:
        """
        :param model: Model to be run across all 'bounds' for 'year'.
        :param year: The year to run the model for.
        :param bounds: A nested list containing lon, lat vertices of boundaries for model to run at. Running over too
            large of a domain can exceed Earth Engine's memory limit so iterating over smaller bounds ensures this won't
            happen.
        :param name: The name of the final asset.
        :param asset_path: The path to the asset directory or image collection to write output to
            (e.g. 'users/username/write_to_this_directory')
        :param kwargs: Dict of additional arguments necessary to run the model.
        """
        self.year = year
        self.bounds = bounds
        self.model = model
        self.kwargs = kwargs
        self.unq_id = unq_id
        self.asset_path = asset_path
        self.scale = scale
        self.out_name = os.path.join(asset_path, self.unq_id+'-'+str(self.year)).replace('\\', '/')
        self.ee_loc = shutil.which('earthengine')
        self.asset_names = self._make_asset_names()

    def _make_asset_names(self) -> List:
        """
        :return: List of names for intermediate mosaic tiles
        """
        check = sp.check_output([self.ee_loc, '--no-use_cloud_api', 'ls', self.asset_path])
        if check.decode() == "No such folder or collection: '{}'.\n".format(self.asset_path):
            sp.call([self.ee_loc, '--no-use_cloud_api', 'create', 'folder', self.asset_path])

        return [os.path.join(x, y).replace('\\', '/') for (x, y) in list(zip([self.asset_path] * len(self.bounds),
                                                                             [self.unq_id+'-intermediate_' + str(z) for
                                                                              z in range(len(self.bounds))]))]

    def _export_tiled_images(self) -> None:
        """
        :return: None, begins tasks for running a model across all tiles
        """
        for i in range(len(self.bounds)):
            print('Starting {} for {}'.format('Intermediate Tile ' + str(i), self.unq_id))
            roi = ee.Geometry.Polygon(self.bounds[i]).bounds()

            dat = self.model(year=self.year, roi=roi, **self.kwargs)
            dat = dat.mean()  # Take the mean across the year
            out = ee.Image(dat).multiply(100).toInt16()  # Convert to int16 to save space

            task = ee.batch.Export.image.toAsset(
                image=out,
                description='{} for {}'.format('Intermediate Tile ' + str(i), self.unq_id),
                scale=self.scale,
                maxPixels=1e13,
                region=self.bounds[i],
                assetId=self.asset_names[i]
            )

            task.start()

    @staticmethod
    def mosaic_tiles(asset_path: str, unq_id: str, scale: int, roi: List = [[[-126.74, 50.06],
                                                                             [-126.74, 22.75],
                                                                             [-65.74, 22.75],
                                                                             [-65.74, 50.06]]]) -> None:
        """
        :param scale:
        :param asset_path: Path to asset imageCollection or folder in earthengine where assets will be written.
        :param unq_id: Unique name identifier for the final mosaiced
        :param roi: Nested list of [lon, lat] that define the bounds of the image to export.
        :return: None, begins EarthEngine task of exporting mosaic to asset.
        """
        check = sp.check_output([shutil.which('earthengine'),
                                 '--no-use_cloud_api', 'ls', asset_path]).decode().split('\n')
        regex = re.compile(unq_id)
        out = ee.ImageCollection.fromImages(
            ee.List([ee.Image(x) for x in list(filter(regex.search, check))])
        ).mosaic()

        print('Mosaicing Intermediate Tiles for {}'.format(unq_id))

        task = ee.batch.Export.image.toAsset(
            image=out,
            description='Mosaicing {}'.format(unq_id),
            scale=scale,
            maxPixels=1e13,
            region=roi,
            assetId=os.path.join(asset_path, unq_id).replace('\\', '/')
        )

        task.start()

    def _mosaic_tiles(self, roi=[[[-126.74, 50.06],
                                  [-126.74, 22.75],
                                  [-65.74, 22.75],
                                  [-65.74, 50.06]]]):
        """
        Class method that implements the static 'mosaic_tiles' method.
        :param roi:
        :return:
        """
        self.mosaic_tiles(self.asset_path, self.unq_id, self.scale, roi)

    def _rm_intermediates(self):
        """
        :return: None, removes all intermediate tiles that were used to create mosaic.
        """
        for asset in self.asset_names:
            sp.call([self.ee_loc, '--no-use_cloud_api', 'rm', asset])

    def run_model(self, wait_time: int, roi: List = [[[-126.74, 50.06],
                                                      [-126.74, 22.75],
                                                      [-65.74, 22.75],
                                                      [-65.74, 50.06]]]):
        """
        Runs all tiles across domain, mosaics them into one continuous image, then removes intermediate tiles
        :param wait_time: Time in seconds to wait between checking whether or not any tasks are running.
        :param roi: Nested list of 4 [lon, lat] locations representing the bounds of the region to run the model. The
            default is the CONUS domain.
        :return:
        """
        self._export_tiled_images()
        while task_running(self.ee_loc):
            print('Cannot mosaic tiles because some intermediate tasks are still running.'
                  ' Waiting {} seconds and will try again.'.format(wait_time))
            time.sleep(wait_time)

        self._mosaic_tiles(roi)
        while task_running(self.ee_loc):
            print('Mosaicing tiles... Waiting {} seconds and will try '
                  'removing intermediate tiles again.'.format(wait_time))
            time.sleep(wait_time)

        self._rm_intermediates()
