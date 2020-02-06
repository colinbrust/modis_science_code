import spotpy
import pandas as pd
import numpy as np
from typing import Dict, Callable, List


class Calibration(object):
    """
    Setup class for model calibration using the spotpy package. Everything implemented in this class is heavily borrowed
    from the Spotpy documentation (https://github.com/thouska/spotpy)..
    """

    def __init__(self, param_dict: Dict, model: Callable, dataset: pd.DataFrame) -> None:

        self.param_dict = param_dict
        self.params = self.gen_params()
        self.model = model
        self.dataset = dataset
        self.param_names = [x for x in self.param_dict]

    def gen_params(self) -> List[spotpy.parameter.Uniform]:
        """
        :return: List of Spotpy Uniform parameters used for model calibration
        """
        return [spotpy.parameter.Uniform(x,
                                         low=self.param_dict[x]['min'],
                                         high=self.param_dict[x]['max'],
                                         optguess=self.param_dict[x]['guess']) for x in self.param_dict]

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector) -> np.array:
        """
        Run model with new random parameter values
        :param vector: List of parameter values returned from gen_params
        :return: np.array containing the new model evaluation generated using new parameters
        """
        param_inputs = [vector[x] for x in self.param_names]
        param_inputs = dict(zip(self.param_names, param_inputs))
        dat = self.model(self.dataset, **param_inputs)
        return dat.simulation.values

    def evaluation(self) -> np.array:
        """
        Get the ground truth values for comparison with model
        :return: np.array containing values that correspond to model simulations
        """
        return self.dataset.target.values

    def objectivefunction(self, simulation: np.array, evaluation: np.array) -> float:
        """
        Compute RMSE error between model and observations
        :param simulation: np.array of model simulations
        :param evaluation: np.array of observations
        :return: RMSE between simulation and evaluation
        """
        objectivefunction = -spotpy.objectivefunctions.rmse(evaluation, simulation)
        return objectivefunction
