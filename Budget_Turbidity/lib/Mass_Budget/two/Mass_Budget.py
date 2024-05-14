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
from diagnostic_function import*
import fluidfoam
#from fluidfoam import OpenFoamSimu
#from Filter_Snapshot import*
#from Budget_Ingradient import*
import sys

#sys.path.append('Operations/read_write')

sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0]))

path = sys.path[0]

pathsrc = os.path.split(path)[0]

sys.path.append(pathsrc)

print ("Manohar in Mass budget = \t", sys.path)

from Operations.spatial_operation.two.Budget_Ingradient import Budget_Integrand
from Operations.read_write.diagnostic_function import write_hdf5
from Operations.extract_time.Filter_Snapshot import Filter_Snapshot

#from diagnostic_function import*

class Mass_Budget:

    def __init__(self, path, dim, budget, parameter, value):
        self.path = path
        self.dim = dim
        self.budget = budget
        self.parameter = parameter
        self.value = value

        print ("Manohar in Mass flux = \t", self.path)
        sol = self.path
    
        BF = Filter_Snapshot(sol, self.dim)

        T       = BF.tprocess
        xbed    = BF.xbed
        Yb      = BF.Yb
        tread   = BF.tread
        Xb = BF.Xb

        n1d = len(xbed)
        Nt = len(T)

        print ("n1d = \t", n1d)
        print ("Nt = \t", Nt)

        time        = np.zeros(Nt)
        dt          = np.zeros(Nt)
        dybeddt_l   = np.zeros((n1d, Nt))
        dybeddt_h   = np.zeros((n1d, Nt))
        balance     = np.zeros((n1d, Nt))
        xf          = np.zeros(Nt)

        xmax = np.amax(xbed)
        ymax = np.amax(Yb[0, :, 0])

        k = -1

        prec = 9

        print ("Path in Mass budget = \t", path)

        f = open('index_zero_xf.d', 'w')

        direc = './output_%s_%s/%s/%s'%(self.dim, self.budget, self.parameter, self.value) 
        
        try:
            os.makedirs(direc, exist_ok = True)
            print("Directory '%s' created successfully" % direc)
        except OSError as error:
            print("Directory '%s' can not be created" % direc) 

        for t in T:

            print ("reading time: %s s"%t)
            k = k + 1
            time[k] = float(t)

            alpha =  readfield(sol, t, 'alpha.a', structured = True, precision = prec )
            Ua    =  fluidfoam.readvector(sol, t, 'U.a', structured = True, precision = prec)

            if t != tread[-1]:
                dt[k] = float(tread[k+1]) - time[k]
                alphadt = readfield(sol, tread[k+1], 'alpha.a', structured = True, precision = prec)
                Uadt      = readfield(sol, tread[k+1], 'U.a', structured = True, precision = prec)

            else:
                dt[k] = dt[k-1]


            xf[k] = np.max(Xb[np.where(alpha>1e-3)])
            range_inte = np.where(alpha[:, 0, 0]>1e-3)[0]
            inte = range_inte[-1]

            range_inte_dt = np.where(alpha[:, 0, 0]>1e-3)[0]
            inte_dt = range_inte_dt[-1]

            xzero = np.where(xbed>0)[0]
            index_xzero =  xzero[0]

            print ("print xzero = \t", index_xzero)
            print ("print inte = \t", inte)

            f.write('%d \t %g \n'%(time[k], xf[k]))

            BI   = Budget_Integrand(xbed, Yb, alpha=alpha)
            BIdt = Budget_Integrand(xbed, Yb, alpha=alphadt)

            if xbed[index_xzero]>=0 and xbed[inte]<=xmax - 0.5*ymax:
                xf[k] = np.max(Xb[np.where(alpha>1e-3)])

                ############### find the bed height

                ybed_l   = np.zeros((n1d))
                ybeddt_l = np.zeros((n1d))
                ybed_h   = np.zeros((n1d))
                ybeddt_h = np.zeros((n1d))

                bed_l   = np.zeros((n1d), dtype = int)
                beddt_l = np.zeros((n1d), dtype = int)
                bed_h   = np.zeros((n1d), dtype = int)
                beddt_h = np.zeros((n1d), dtype = int)

                bed_l, bed_h, ybed_l, ybed_h         = BI.bed_height()
                beddt_l, beddt_h, ybeddt_l, ybeddt_h = BIdt.bed_height()

                ########### derivative with space
                
                dybdx_l = BI.x_deri(ybed_l, xbed)
                dybdx_h = BI.x_deri(ybed_h, xbed)

                dybdxdt_l = BIdt.x_deri(ybeddt_l, xbed)
                dybdxdt_h = BIdt.x_deri(ybeddt_h, xbed)

                ########## derivative with time

                dybdt_l = np.zeros(n1d)
                dybdt_h = np.zeros(n1d)

                dybdt_l = (ybeddt_l - ybed_l) / dt[k]
                dybdt_h = (ybeddt_h - ybed_h) / dt[k]

                dybeddt_l[:, k] = BI.SurfacevalueMul(alpha[:, :, :], dybdt_l, bed_l)
                dybeddt_h[:, k] = BI.SurfacevalueMul(alpha[:, :, :], dybdt_h, bed_h)

                ######### nonlinear terms

                dybdx_uphi_l = np.zeros(n1d)
                dybdx_uphi_h = np.zeros(n1d)

                vphi_l = np.zeros(n1d)
                vphi_h = np.zeros(n1d)

                uphi = alpha*Ua[0, :, :, :]
                vphi = alpha*Ua[1, :, :, :]

                uphidt = alphadt*Uadt[0, :, :, :]
                vphidt = alphadt*Uadt[1, :, :, :]
                
                dybdx_uphi_l = (BI.SurfacevalueMul(uphi, dybdx_l, bed_l) + BIdt.SurfacevalueMul(uphidt, dybdxdt_l, beddt_l)) / 2.0
                dybdx_uphi_h = (BI.SurfacevalueMul(uphi, dybdx_h, bed_h) + BIdt.SurfacevalueMul(uphidt, dybdxdt_h, beddt_h)) / 2.0

                vphi_l = (BI.Surfacevalue(vphi, bed_l) + BIdt.Surfacevalue(vphidt, beddt_l)) / 2.0
                vphi_h = (BI.Surfacevalue(vphi, bed_h) + BIdt.Surfacevalue(vphidt, beddt_h)) / 2.0

                ############## Storage

                Inte_phi   = np.zeros(n1d)
                Inte_phidt = np.zeros(n1d)
                storage    = np.zeros(n1d)

                Inte_phi   = BI.Integrate_flux(alpha, Yb, bed_l, bed_h)
                Inte_phidt = BIdt.Integrate_flux(alphadt, Yb, beddt_l, beddt_h)
                Storage    = (Inte_phidt - Inte_phi) / dt[k]

                ############# flux

                Inte_uphi = np.zeros(n1d)
                uphiflux  = np.zeros(n1d)
                Inte_phi = BI.Integrate_flux(uphi, Yb, bed_l, bed_h) + BIdt.Integrate_flux(uphidt, Yb, beddt_l, beddt_h)
                uphiflux  = BI.x_deri(Inte_phi, xbed)

                ############# Writing files

                write_hdf5( direc + '/' + 'Storage_%g.h5'% (time[k]), Storage)
                write_hdf5( direc + '/' + 'Inte_phi_%g.h5'% (time[k]), Inte_phi)
                write_hdf5( direc + '/' + 'Inte_phidt_%g.h5'% (time[k]), Inte_phidt)

                ########### flux
                write_hdf5(direc + '/' + 'dybdx_uphi_l_%g.h5'% (time[k]), dybdx_uphi_l)
                write_hdf5(direc + '/' + 'dybdx_uphi_h_%g.h5'% (time[k]), dybdx_uphi_h)
                write_hdf5(direc + '/' + 'vphi_l_%g.h5'%(time[k]), vphi_l)
                write_hdf5(direc + '/' + 'vphi_h_%g.h5'% (time[k]), vphi_h)
                write_hdf5(direc + '/' + 'uphiflux_%g.h5'% (time[k]), uphiflux)

                ############# bottom part

            else:
                break

        write_hdf5(direc + '/' + 'xbed.h5', xbed)
        write_hdf5(direc + '/' + 'dybeddt_l.h5', dybeddt_l)
        write_hdf5(direc + '/' + 'dybeddt_h.h5', dybeddt_h)


