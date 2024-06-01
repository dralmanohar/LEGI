# import python packages
from fluidfoam import readmesh,readfield
import os
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import CubicSpline
import subprocess
import sys
import matplotlib.pyplot as plt
from netCDF4 import Dataset
#from diagnostic_function import*
import fluidfoam
#from fluidfoam import OpenFoamSimu
#from Operations/extract_time.Operations/extract_time import*
import sys

sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))

path = sys.path[0]

pathsrc = os.path.split(path)[0]

sys.path.append(pathsrc)

print ("Path inside Momentum Budget = \t", sys.path)

from Operations.spatial_operation.two.Budget_Ingradient import Budget_Integrand
from Operations.read_write.diagnostic_function import write_hdf5
from Operations.extract_time.Filter_Snapshot import Filter_Snapshot

class Momentum_Budget:
	def __init__(self, flow_type, path, dim, budget, parameter, value):
		
		self.flow_type = flow_type
		self.path = path
		self.dim = dim
		self.budget = budget
		self.parameter = parameter
		self.value = value
		if self.flow_type == 'Poiseuille':
			from two.Poiseuille import Momentum_Budget_Poiseuille 
			Momentum_Budget_Poiseuille(self.path, self.dim, self.budget, self.parameter, self.value)
		elif self.flow_type == 'Turbidity':
			from Turbidity import Momentum_Budget_Turbidity
			Momentum_Budget_Turbidity(self.path, self.dim, self.budget, self.parameter, self.value)
        

