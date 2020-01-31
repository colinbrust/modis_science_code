import pandas as pd
import spotpy
from scripts.Dataset import Dataset
from typing import Dict, Callable


class Calibration(object):

    def __init__(self, param_dict: Dict, model: Callable, dataset: Dataset) -> None:

        self.param_dict = param_dict
        self.params = self.gen_params()
        self.model = model
        self.dataset = dataset
        self.param_names = [x for x in self.param_dict]

    def gen_params(self):

        return [spotpy.parameter.Uniform(x,
                                         low=self.param_dict[x]['min'],
                                         high=self.param_dict[x]['max'],
                                         optguess=self.param_dict[x]['guess']) for x in self.param_dict]

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):

        param_inputs = [vector[x] for x in self.param_names]
        param_inputs = dict(zip(self.param_names, param_inputs))
        simulations = self.model(self.dataset.df, **param_inputs)

        return simulations.simulated_values

    def evaluation(self):
        return self.dataset.target_values

    def objectivefunction(self, simulation, evaluation):
        objectivefunction = -spotpy.objectivefunctions.rmse(evaluation, simulation)
        return objectivefunction
