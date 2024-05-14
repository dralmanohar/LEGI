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
#from Budget_function import *

#path = '/.fsnet/data/legi/calcul9/home/sharma6ma/useful/project/17LES_SHEET/manohar/palagram/run_37_theta_0/Conference_EGU/2D/250_um/theta_0'
#path = '/.fsnet/data/legi/calcul9/home/sharma6ma/useful/project/17LES_SHEET/manohar/Channel/LES/Poisuelli_flow' #run_37_theta_0/Conference_EGU/2D/250_um/theta_0
path = '/.fsnet/data/legi/calcul9/home/sharma6ma/useful/project/17LES_SHEET/manohar/palagram/run_37_theta_0/Conference_EGU/2D/250_um/theta_0_pressure_outlet_0'

#path = '/.fsdyn_people/sharma6ma/project/21PALAGRAM/Manohar/3D_Poiseuilli_flow'
file_name = "configurationsMarie_run07.csv"

if __name__ == '__main__':

    ######### define parse

    parser = argparse.ArgumentParser()
    
    ########### Add element 
    parser.add_argument("Budget", help = 'Give the name of budget', choices = ["Movie", "Mass", "Momentum"])
    parser.add_argument("dim", help = "Give dimension")
    parser.add_argument("parameter", help = 'Give the parameter', choices = ["diameter", "angle", "phi"])
    parser.add_argument("value", help = 'Give the value of the parameter')
    #parser.add_argument("dim", help = 'Give dimension')
    
    ########## finalize
    args = parser.parse_args()
    
    if args.Budget=="Momentum":
    
        print ("Manohar in compute file")
        sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))
        pathsrc = sys.path[0]
        path_lib = os.path.join(pathsrc, 'lib')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
        sys.path.append(path_lib)

        print ("sys path Momentum = \t", sys.path)

        from Budget_function import Budget_seperator

        Momentum = Budget_seperator(path, args.Budget, args.dim, args.parameter, args.value)
        #Momentum = Momentum_Budget(path, args.Budget, args.dim, args.parameter, args.value)
        print ("Manohar in Momentum = \t", args.Budget, args.dim, args.parameter, args.value)
    
    elif args.Budget=="Mass":
        
        sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))
        pathsrc = sys.path[0]
        path_lib = os.path.join(pathsrc, 'lib')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
        sys.path.append(path_lib)
        
        from Budget_function import Budget_seperator
        #print ("Mass budget in compute= \t", path)
        Mass = Budget_seperator(path, args.Budget, args.dim, args.parameter, args.value)
        print ("Manohar in Mass flux = \t", args.Budget, args.dim, args.parameter, args.value)
    elif args.Budget=="Movie":
        
        sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))
        pathsrc = sys.path[0]
        path_lib = os.path.join(pathsrc, 'lib')#sys.path[0] + '/' + 'Momentun_Budget' #+ '/' + 'two'
        sys.path.append(path_lib)
        
        from Budget_function import Budget_seperator
        
        #print ("Manohar in Make Movie = \t", args.Budget, args.dim, args.parameter, args.value)
        Movie = Budget_seperator(path, args.Budget, args.dim, args.parameter, args.value, file_name = file_name)
        #Movie = Make_Movie(path, file_name, args.dim)
        #print ("Manohar in Make Movie = \t", args.Budget, args.dim, args.parameter, args.value)




