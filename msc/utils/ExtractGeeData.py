import ee
import os
import pandas as pd
from typing import Callable, TypeVar, Generic

# Area in meters to buffer around each point and run the model at
BUFFER_SIZE = 10000
T = TypeVar('T')


class ExtractGeeData(object):
    """
    Class to extract GEE data at flux tower locations in order to run and calibrate MOD16 and MOD17.
    """
    def __init__(self, template: str, model: Callable, out_dir: str, **kwargs: Generic[T]) -> None:
        """
        :param template: File path to a .csv containing columns 'name', 'date', 'lat', 'longitude', and 'target', which
        correspond to the site name where ground observations came from, the date of the observation, the latitude and
        longitude of the site, and the recorded observation from the site, respectively.

        :param model: A function that runs the desired model at the site locations specified in the template.

        :param out_dir: Directory where data will be written to

        :param kwargs: Keyword arguments for the model that can't be derived from the template file. For example, to run
        MOD16, daylength, elevation, and meteorology imageCollections would be necessary, as the roi and year arguments
        can be taken from the template file.
        """
        self.template = pd.read_csv(template)
        self.model = model
        self.kwargs = kwargs
        self.out_dir = out_dir
        self.full_dataset = pd.DataFrame()
        self.run_instructions = self._get_locations_years()
        self.restart_instructions = pd.DataFrame()

    def _get_locations_years(self) -> pd.DataFrame:

        """
        :return: Pandas DataFrame containing a unique combination of all sites and years that will be used to direct
        where and when to run the model
        """
        df = self.template.copy()

        df['year'] = df.date.str[0:4]
        df['year'] = df['year'].astype(int)
        df = df[['name', 'year', 'lat', 'lon']]
        df = df.drop_duplicates()

        return df

    @staticmethod
    def _reducer(img: ee.Image, point: ee.geometry.Geometry, name: str) -> ee.Feature:
        """
        Reduces an area around a point to a single value and returns as an ee.Feature
        :param img: ee.Image to reduce
        :param point: ee.geometry.Geometry that represents the point to reduce in img
        :param name: The site name of 'point'
        :return: ee.Feature containing the reduced value
        """
        reduced = img.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=point,
            scale=500,
            maxPixels=1e8
        )

        time = img.get('system:time_start')

        return ee.Feature(None, reduced).set('system:time_start', time).set('name', name)

    def _run_single_location(self, name: str, year: int, lat: float, lon: float) -> pd.DataFrame:
        """
        :param name: Name of the site where the model will be run
        :param year: Year that the model will be run for
        :param lat: Latitude of site location
        :param lon: Longitude of site location
        :return: Pandas data frame containing model results
        """
        pnt = ee.Geometry.Point([lon, lat]).buffer(BUFFER_SIZE)
        model_results = self.model(roi=pnt, year=year, **self.kwargs)
        model_results = model_results.map(lambda img: self._reducer(img=img, point=pnt, name=name))
        model_results = model_results.getInfo()

        listed = model_results['features']

        df = []
        for day in listed:
            df.append(day['properties'])

        return pd.DataFrame(df)

    def _restart(self) -> None:
        """
        Whether or not the model run should be restarted using temp_save.csv
        :return: None
        """
        completed = pd.read_csv(os.path.join(self.out_dir, 'temp_save.csv'))
        completed = completed[['name', 'year']].drop_duplicates()

        for _, row in completed.iterrows():
            self.run_instructions = self.run_instructions[(self.run_instructions['year'] != row['year']) |
                                                          (self.run_instructions['name'] != row['name'])]

    def run_model(self, restart: bool = False) -> pd.DataFrame:
        """
        :param restart: If True, will look for temp_save.csv file in data_dir and update the list of sites/years at
        which to run the model according to what was already done in temp_save.csv
        :return: Pandas DataFrame containing model results for all sites and locations. Also saves df out as csv to
        out_dir
        """
        if restart:
            self._restart()

        for _, row in self.run_instructions.iterrows():

            print('Running model at {} for {}'.format(row['name'], row['year']))

            # Because the model is pretty computationally demanding, it can exceed EarthEngine's memory limit.
            # This prevents the code from crashing if that happens.
            while True:
                try:
                    df = self._run_single_location(name=row['name'], year=row['year'],
                                                   lat=row['lat'], lon=row['lon'])

                except ee.ee_exception.EEException as e:
                    print(e)
                    print('Re-running model at {} for {}'.format(row['name'], row['year']))
                    continue
                break

            df['year'] = row['year']
            self.full_dataset = pd.concat([self.full_dataset, df])

            # Save out data every model run in case of crash.
            self.full_dataset.to_csv(os.path.join(self.out_dir, 'temp_save.csv'), index=False)

        self.full_dataset.to_csv(os.path.join(self.out_dir, 'modeled_results.csv'), index=False)
        return self.full_dataset
