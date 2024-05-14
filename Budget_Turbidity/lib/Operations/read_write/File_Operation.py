# Import python packages
import subprocess,sys,os
import csv
import fluidfoam 

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import scipy.optimize as so
import scipy.stats as ss
from netCDF4 import Dataset
#from diagnostic_function import*



class File_Operation:
    def __init__(self, csvfile):
        self.csvfile = csvfile

    def decomment(self):
        for row in self.csvfile:
            self.raw = row.split('#')[0].strip()
            if self.raw: yield self.raw

    def create_variable(self, netcdf_group, name, data, dimensions=None, std=None, unit=None, comments=None, type='float64'):
        if dimensions is not None:
            self.var = netcdf_group.createVariable(name, type, (dimensions))
        else:
            self.var = netcdf_group.createVariable(name, type)

        self.var[:] = data
        if std is not None:
            self.var.std = std
        if unit is not None:
            self.var.unit = unit
        if comments is not None:
            self.var.comments = comments

    def fix_unit(self):
        if self.var.unit == 'cm':
            factor = 1e-2
            newunit = 'm'
        elif var.unit == 'cm2':
            factor = 1e-4
            newunit = 'm2'
        elif self.var.unit == 'g/cm3':
            factor = 1e3
            newunit = 'kg/m3'
        elif self.var.unit == 'cm/s':
            factor = 1e-2
            newunit = 'm/s'
        elif self.var.unit == '[cm, cm ,cm]':
            factor = 1e2
            newunit = '[m, m ,m]'
        elif self.var.unit == 'mum':
            factor = 1e-6
            newunit = 'm'
        if 'factor' in locals():
            self.var[:] = self.var[:]*factor
            self.var.unit = newunit
            if hasattr(var, 'std'):
                self.var.std = self.var.std*factor





