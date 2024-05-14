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
#from Filter_Snapshot import*
#from Budget_Ingradient import*
import argparse
#from Mass_Budget import*
#from Momentum_Budget import*
#from Make_Movie import*
import sys

class Budget_seperator:
    def __init__(self, path,  budget, dim, parameter, value, file_name = None):
        
        self.path = path
        self.dim = dim
        self.budget = budget
        self.dim = dim
        self.parameter = parameter
        self.value = value
        self.file_name = file_name
        
        if self.budget=='Mass':
            if self.dim == '2D':
                
                sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))
                pathsrc = sys.path[0]

                #print ("pathsrc for the file = \t", pathsrc)
                path_libL = os.path.join(pathsrc, 'lib')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
                path_libS = os.path.join(path_libL, 'Mass_Budget')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
                path_lib = os.path.join(path_libS, 'two') 

                sys.path.append(path_lib)
                
               # sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))
               # pathsrc = sys.path[0]
               # path_lib = os.path.join(pathsrc, 'Mass_Budget')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
               # path_lib = os.path.join(path_lib, 'two') 

               # sys.path.append(path_lib)
                
                from Mass_Budget import Mass_Budget

                self.Mass = Mass_Budget(self.path, self.dim, self.budget, self.parameter, self.value)
            
            elif self.dim == '3D':
                sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))
                pathsrc = sys.path[0]

                #print ("pathsrc for the file = \t", pathsrc)
                path_libL = os.path.join(pathsrc, 'lib')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
                path_libS = os.path.join(path_libL, 'Mass_Budget')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
                path_lib = os.path.join(path_libS, 'three') 

                sys.path.append(path_lib)
                
                #sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))
                #pathsrc = sys.path[0]
                #path_lib = os.path.join(pathsrc, 'Mass_Budget')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
                #path_lib = os.path.join(path_lib, 'three') 

                #sys.path.append(path_lib)
                
                from Mass_Budget import Mass_Budget

                self.Mass = Mass_Budget(self.path, self.dim, self.budget, self.parameter, self.value)

#                sys.path.append('./Mass_Budget/')
 #               from three.Mass_Budget import Mass_Budget
  #              self.Mass = Mass_Budget(path, budget, parameter, value)
        elif self.budget == 'Momentum':
            
            if self.dim == '2D':
                sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))
                pathsrc = sys.path[0]

                #print ("pathsrc for the file = \t", pathsrc)
                path_libL = os.path.join(pathsrc, 'lib')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
                path_libS = os.path.join(path_libL, 'Momentum_Budget')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
                path_lib = os.path.join(path_libS, 'two') 

                sys.path.append(path_lib)
                
                print ("path for the 2D momentum budget = \t", sys.path)

                from Momentum_Budget import Momentum_Budget
                self.Momentum = Momentum_Budget(self.path, self.dim, self.budget, self.parameter, self.value)
            
            elif self.dim == '3D':
                sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))
                pathsrc = sys.path[0]

                #print ("pathsrc for the file = \t", pathsrc)
                path_libL = os.path.join(pathsrc, 'lib')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
                path_libS = os.path.join(path_libL, 'Momentum_Budget')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
                path_lib = os.path.join(path_libS, 'three') 

                sys.path.append(path_lib)

    #            sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))
#                pathsrc = sys.path[0]
 #               path_lib = os.path.join(pathsrc, 'Momentum_Budget')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
  #              path_lib = os.path.join(path_lib, 'three')
   #             sys.path.append(path_lib)
                
                from Momentum_Budget import Momentum_Budget
                
                self.Momentum = Momentum_Budget(self.path, self.dim, self.budget, self.parameter, self.value)
        
        elif self.budget == 'Movie':
                sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))
                pathsrc = sys.path[0]

                #print ("pathsrc for the file = \t", pathsrc)
                path_libL = os.path.join(pathsrc, 'lib')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
                path_lib = os.path.join(path_libL, 'Front_Movie')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'

                sys.path.append(path_lib)

                #sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))
                #pathsrc = sys.path[0]
                #path_lib = os.path.join(pathsrc, 'Front_Movie')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
                #sys.path.append(path_lib)
                
                #print ("path for the movie library = \t", sys.path)
                from Front_Movie import Make_Movie
                #print ("filename in Make movie = \t", file_name)
                self.Movie = Make_Movie(self.path, self.file_name)#, self.budget, self.parameter, self.value)

            
#                sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))
#                path = sys.path[0]
#                path_lib = os.path.join(path, 'Make_Movie')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
                
#                print ("Path for Making movie = \t", path_lib)
#                sys.path.append(path_lib)

#                from Make_Movie import Make_Movie

#                if file_name!=None:
#                    print ("Manohar in Budget for Making Movie")
#                    self.Movie = Make_Movie(self.path, file_name = self.file_name)
#                else:
#                    print ("Please provide the name of the csv file which has the cases")
    




