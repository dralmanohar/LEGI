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
import fluidfoam
#from fluidfoam import OpenFoamSimu


class Filter_Snapshot:

    def __init__(self, path, dim):
        self.path = path
        self.dim = dim

        prec = 13

        sol = self.path 
        
        print ("path in filter operation = \t", sol)

        Xb, Yb, Zb = readfield(sol, 'constant', 'C', structured = True, precision = prec)

        (n1d, ny, n2d) = np.shape(Xb)
    
        self.ny  = ny
        #self.xbed = Xb[:, 0, 0]
        self.Yb = Yb
        self.Xb = Xb
        #zbed = Zb[0, 0, :]

        if dim =='2D':
            self.xbed = Xb[:, 0, 0]
            self.zbed = Zb[0, 0, :]
        elif dim == '3D':
            self.xbed = Xb[:, 0, :]
            self.zbed = Zb[:, 0, :]
            print ("Shape in Filter snapshot xbed = \t", Xb[:, 0, :].shape)
        
        try:
            proc = subprocess.Popen(['foamListTimes', '-case', sol], stdout = subprocess.PIPE)

        except:
            print ("foamListTime: command not found")
            print ("Do you have loaded OpenFoam environment")
            sys.exit(0)

        output = proc.stdout.read()

        tread = output.decode().rstrip().split('\n')

        tread.insert(0, '0')
        del tread[0]

        self.tprocess = tread.copy()
        self.tread = tread
        del self.tprocess[-1]

