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
from Compute_Front import*

sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))

path = sys.path[0]

pathsrc = os.path.split(path)[0]

sys.path.append(pathsrc)

#print ("Manohar in Make movie = \t", sys.path)

#from Make_Movie.Compute_Front import Compute_Front
from Operations.spatial_operation.two.Budget_Ingradient import Budget_Integrand
from Operations.spatial_operation.three.Budget_Ingradient import Budget_Integrand
from Operations.read_write.diagnostic_function import write_hdf5
from Operations.read_write.File_Operation import File_Operation
from Operations.extract_time.Filter_Snapshot import Filter_Snapshot



class Make_Movie:
    def __init__(self, path, file_name, dim = None):
        self.path = path
        self.file_name = file_name
        self.dim = dim



        ########################################################################
        # Config name
        ########################################################################

        config_name = self.path         #'./'
        config_file = self.file_name    #'configurationsMarie_run07.csv'

        makeMovies = True

 #       print ("Manohar in Make movie class")


        # Get the current working directory
        pyPath = self.path #os.getcwd()
        
        print ("basepath before = \t", pyPath)
        basePath = os.path.split(pyPath)[0]
        #basePath = os.path.split(basePath)[0]

        print ("base path after = \t",basePath)

        #print("basePath = \t", basePath)

        # Read csv file containing simulations parameters
        name_simu = [ ]
        caseDictList = [ ]

        #FO = file_operation()

        with open(config_file, newline='') as csvfile:
            FO = File_Operation(csvfile)
            configurations = csv.reader(FO.decomment())
            #print ("configurations = \t",configurations)
            for line in configurations:
                #print ("line = \t", line)
                row=line[0].split(';')
                #print ("row = \t", row)
                caseName = row[0].strip()
                name_simu.append(caseName)
                case = {
			        "runName": caseName,
			        "phi" : float(row[1]),
			        "alpha" : float(row[2]),
			        "gsinalpha": float(row[3]),
			        "gcosalpha": float(row[4]),
		    	    "rho_f" : float(row[5]),
			        "rho_p": float(row[6]),
			        "L0" : float(row[7]),
			        "H0" : float(row[8]),
			        "W0" : float(row[9]),
			        "partDiameter": float(row[10]),
			        "v_s": float(row[11]),
		            }
                caseDictList.append(case)
    
        basePathoutput = os.path.join(pyPath,config_name)
        
        pyPathout = os.getcwd()
        print ("Manohar Make movie output = \t", pyPathout)
        basePathoutput = pyPathout#os.path.split(pyPathout)
        outputPath = os.path.join(basePathoutput,'output')

        print ("outputPath Make movie = \t", outputPath)
        if os.path.exists(outputPath)==False:
            os.mkdir(outputPath)
        
        for case in caseDictList:
            print ("case = \t", case)
            caseName = case["runName"]
            print('post-process : ', caseName)
            compute_front = Compute_Front(basePath,config_file,outputPath,case,makeMovies, dim = '2D')


