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
#from Filter_Snapshot import*
#from Budget_Ingradient import*
#from file_operation import*


sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))

path = sys.path[0]

pathsrc = os.path.split(path)[0]

sys.path.append(pathsrc)

#print ("Manohar in Mass budget = \t", sys.path)

from Operations.spatial_operation.two.Budget_Ingradient import Budget_Integrand
from Operations.spatial_operation.three.Budget_Ingradient import Budget_Integrand
from Operations.read_write.diagnostic_function import write_hdf5
from Operations.read_write.File_Operation import File_Operation
from Operations.extract_time.Filter_Snapshot import Filter_Snapshot



class Compute_Front:

    def __init__(self,basePath, config_file, outputPath, case, makeMovies, dim=None ):
        self.basePath = basePath
        self.config_file = config_file
        self.outputPath = outputPath
        self.case = case
        self.makeMovies = makeMovies
        self.dim = dim
        #self.csvfile = csvfile
        
        caseName =self. case["runName"].split("/")[-1]
        
        print ("caseName for compute front = \t", caseName)

        sol = os.path.join(self.basePath, caseName)
    
        path = sol+ '/' + 'images'

        print ("sol for Manohar = \t", sol)
        print ("Path for images for Manohar = \t", path)


        try:
            os.mkdir(path)
            print("Folder %s created!" % path)
        except FileExistsError:
            print("Folder %s already exists" % path)

        figPath = os.path.join(sol,'images')

        print ("figPath = \t", figPath)

        figName = 'img'
        frontFile = caseName.split('/')[-1] +'_front'

        print ("fron file in the function = \t", frontFile)
        rad2deg = 180./np.pi

        nlevels = 100

        levels = np.linspace(0,0.02,nlevels)
        cbar_levels = np.linspace(0,0.02,6)
        cbar_labels = ['0','0.004','0.008','0.012','0.016','0.02']

        print ("caseName = \t",caseName)

        frontfile_name = frontFile.split('/')

        fileName = os.path.join(self.outputPath,caseName)

        print ("save file:",fileName)

        try:
            os.makedirs(fileName)
            print("Folder %s created!" % fileName)
        except FileExistsError:
            print("Folder %s already exists" % fileName)
   
        #if os.path.isfile(fileName):
        #   proc = subprocess.Popen(['rm', '-rf', fileName ], stdout=subprocess.PIPE)

        try:
            proc = subprocess.Popen(['foamListTimes', '-case', sol, '-withZero'], stdout=subprocess.PIPE)
        except:
            print("foamListTimes : command not found")
            print("Do you have load OpenFoam environement?")
            sys.exit(0)

        output = proc.stdout.read()
        tread = output.decode().rstrip().split('\n')
        #print ("Manohar tread = \t", tread)
 

        del tread[-1]
        Nt = len(tread)
        time = np.zeros(Nt)

        print ("Nt = \t", Nt)

        if self.dim=='3D':
            X, Y, Z = fluidfoam.readmesh(sol,structured=True,precision = 9)
            Nx, Ny, Nz = np.shape(X)
        else:
            print ("Manohar in 2D loop")
            X, Y, Z  = fluidfoam.readmesh(sol,structured=True,precision = 6)#9)
            Nx, Ny, Nz  = np.shape(X)
        xmax = np.max(X)
        ymax = np.max(Y)

        try:
            proc = subprocess.Popen(['foamListTimes', '-case', sol, '-withZero'], stdout=subprocess.PIPE)
        except:
            print("foamListTimes : command not found")
            print("Do you have load OpenFoam environement?")
            sys.exit(0)
            print ("Manohar in the compute front loop")
            output = proc.stdout.read()
            tread = output.decode().split('\n')[:-1]
#           print ("Manohar after tread 1 = \t",tread)


        nt=np.size(tread)
        xf=np.zeros(nt)
        yf=np.zeros(nt)
        time=np.zeros(nt)
        print ("nt = \t", nt)

        k = 0
        range_index = 10

#       print ("tread Manohar = \t", tread)

        while k < nt:
            time[k]=float(tread[k])
            #print ("Time Manohar = \t", time[k])
            timename=tread[k]
            alpha_a = fluidfoam.readscalar(sol, timename, 'alpha.a', True,precision=6)
            Ub = fluidfoam.readvector(sol, timename, 'U.b', True,precision=6)
            #Ua = fluidfoam.readscalar(sol, timename, 'alpha.a', True,precision=6)
            try:
                xf[k]=np.max(X[np.where(alpha_a>1e-3)])

            except ValueError:
                pass

            if makeMovies == True:
                print ("time in the ploting loop = \t", time[k], "\t k = \t",k)
                fig = plt.figure(num=k, figsize=(12, 4), dpi=100, facecolor='w', edgecolor='w')
                plt.title('t=' + str(time[k]) + ' s')
                alpha_a[np.where(alpha_a<1e-6)] = 1e-6

                #CS_colors = plt.contourf(X[:,:,0], Y[:,:,0], np.log10(alpha_a[:,:,0]),
                #levels, cmap='gray')
            
                if dim=='3D':
                    CS_colors = plt.contourf(X[:,:,Nz//2], Y[:,:,Nz//2], alpha_a[:,:,Nz//2],levels, cmap='gray')
                else:
                    CS_colors = plt.contourf(X[:,:, 0], Y[:,:, 0], alpha_a[:,:, 0],levels, cmap='gray')
                cb = fig.colorbar(CS_colors)
                cb.set_ticks(cbar_levels)
                cb.set_ticklabels(cbar_labels)
                cb.set_label(r'$\phi$')
                front = plt.plot([xf[k], xf[k]], [0, np.max(Y)],'-r',lw=2)
                plt.axis([np.min(X), np.max(X),np.min(Y), np.max(Y)])
                plt.ylabel('y (m)')
                plt.xlabel('x (m)')

                if len(str(k))>=4:
                    kName = str(k)
                elif len(str(k))>=3:
                    kName = '0'+str(k)
                elif len(str(k))>=2:
                    kName = '00'+str(k)
                else:
                    kName = '000'+str(k)
                plt.savefig(figPath+'/'+figName+str(kName)+'.png',format='png',dpi = 300)
                plt.close(k)
                k+=1

            if self.makeMovies == True:
                #ans1 = input('Do you want to create an .mp4 file ? y=yes / n=no \n')
                ans1 = 'y'
                if ans1=='y':
                    print('Generate the movie using ffmpeg')
                    proc = subprocess.Popen(['rm', '-rf', figPath+'/'+figName+'*.swf'], stdout=subprocess.PIPE)
                    proc = subprocess.Popen(['ffmpeg', '-framerate', '4', '-i', figPath+'/'+figName+'%04d.png', '-c:v','libx264', '-y','-pix_fmt', 'yuv420p', outputPath+'/'+ caseName + '/'+ caseName.split('/')[-1] +'.mp4'], stdout=subprocess.PIPE)
                    (output, err) = proc.communicate()
                    p_status = proc.wait()
        vf = np.zeros(nt)
        dum = np.diff(xf)/np.diff(time)
        vf[:-1] = dum[:]
        vf[-1]  = dum[-1]
        keep_data = np.where(xf<=xmax-0.5*ymax)
        print("keep data = \t",keep_data)

        fileName = os.path.join(fileName,frontFile+'.nc')
   # fileName = os.path.join(fileName,frontFile+'.nc')

        print ("save file:",fileName)

    # if os.path.isfile(fileName):
    #     proc = subprocess.Popen(['rm', '-rf', fileName ], stdout=subprocess.PIPE)

        #NetCFD file creations

        rootgrp = Dataset(fileName, 'w')
        #Dimensions creation
        Ntnc = np.size(keep_data)
        rootgrp.createDimension('time', Ntnc)
        rootgrp.particle_type = 'glass beads'
        #rootgrp.label = run_global
        rootgrp.author = 'Julien'
        rootgrp.lab = 'LEGI'

        variables = [('phi', 'Volume_fraction', '-'),
                ('d', 'Grain size', 'm'),
                ('rho_p', 'Grain density', 'kg/m3'),
                ('v_s', 'settling_velocity', 'm/s'),
                ('alpha', 'slope', 'deg.'),
                ('rho_a', 'Water density', 'kg/m3'),
                ('rho_f', 'Water density', 'kg/m3'),
                ('Current density', 'Current density', 'kg/m3'),
                ]
    
        dico = {}
        FO = File_Operation(self.config_file)

        for (name, key, unit) in variables:
            if name == 'phi':
                dico[key] = case["phi"]
            elif name == 'd':
                dico[key] = case["partDiameter"]
            elif name == 'rho_p':
                dico[key] = case["rho_p"]
            elif name == 'v_s':
                dico[key] = case["v_s"]
            elif name == 'alpha':
                dico[key] = case["alpha"]
            elif name in ['rho_a', 'rho_f']:
                dico[key] = case["rho_f"]
            elif name == 'Current density':
                current_density = dico['Water density']*(1 + dico['Volume_fraction']*(
                    dico['Grain density']-dico['Water density'])/dico['Water density'])
                dico[key] = current_density
            
            FO.create_variable(rootgrp, name,dico[key],std=np.NaN,unit=unit)
        
       # FO = file_operation(self.csvfile)

        FO.create_variable(rootgrp, 'L0', case["L0"],std=np.NaN, unit='m')
        FO.create_variable(rootgrp, 'W0', case["W0"],std=np.NaN, unit='m')
        FO.create_variable(rootgrp, 'H0', case["H0"], std=np.NaN, unit='m')
    
        FO.create_variable(rootgrp, 't', time[keep_data], dimensions=('time'), unit='s')
        FO.create_variable(rootgrp, 'x_front', xf[keep_data], dimensions=('time'), unit='m')
        FO.create_variable(rootgrp, 'v_front', vf[keep_data], dimensions=('time'), unit='m/s')
        rootgrp.close()

