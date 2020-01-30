import ee
import pandas as pd

ee.Initialize()
# Area in meters to buffer around each point and run the model at
BUFFER_SIZE = 500


class ExtractGeeData(object):

    def __init__(self, template, model, **kwargs):
        """
        :param template: File path to a .csv containing columns 'name', 'date', 'lat', 'longitude', and 'target', which
        correspond to the site name where ground observations came from, the date of the observation, the latitude and
        longitude of the site, and the recorded observation from the site, respectively.

        :param model: A function that runs the desired model at the site locations specified in the template.

        :param kwargs: Keyword arguments for the model that can't be derived from the template file. For example, to run
        MOD16, daylength, elevation, and meteorology imageCollections would be necessary, as the roi and year arguments
        can be taken from the template file.
        """
        self.template = pd.read_csv(template)
        self.model = model
        self.full_dataset = pd.DataFrame()
        self.run_instructions = self._get_run_locations_years()

    def _get_run_locations_years(self):

        """
        :return: Pandas DataFrame containing a unique combination of all sites and years that will be used to direct
        where and when to run the model
        """
        df = self.template.copy()

        df['year'] = df.date.str[0:4]
        df = df[['name', 'year', 'lat', 'lon']]
        df = df.drop_duplicates()

        return df

    @staticmethod
    def _reducer(img, point, name):
        reduced = img.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=point,
            scale=500,
            maxPixels=1e8
        )

        time = img.get('system:time_start')

        return ee.Feature(None, reduced).set('system:time_start', time).set('name', name)

    def _run_single_location(self, name, year, lat, lon):

        pnt = ee.Geometry.Point([lon, lat]).buffer(BUFFER_SIZE)
        model_results = self.model(roi=pnt, year=year, **self.kwargs)
        model_results = model_results.map(lambda img: self._reducer(img=img, point=pnt, name=name))
        model_results = model_results.getInfo()