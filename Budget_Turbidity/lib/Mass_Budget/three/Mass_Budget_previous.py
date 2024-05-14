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
from Filter_Snapshot import*
from Budget_Ingradient import*


class Mass_Budget:

    def __init__(self, path, dire, parameter, value):
        self.path = path
        self.dire = dire
        self.parameter = parameter
        self.value = value

    
        sol = self.path
    
        BF = filter_snapshot(sol)

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

        sol = path
        f = open('index_zero_xf.d', 'w')

        direc = './output_2D_%s/%s/%s'%(self.dire, self.parameter, self.value) 
        
        try:
            os.makedirs(direc, exist_ok = True)
            print("Directory '%s' created successfully" % direc)
        except OSError as error:
            print("Directory '%s' can not be created" % direc) 

        for t in T:

            print ("reading time: %s s"%t)
            k = k + 1
            time[k] = float(t)

            alpha =  readfield(sol, t, 'alpha.a', structured = True, precision = 9 )
            Ua    =  fluidfoam.readvector(sol, t, 'U.a', structured = True, precision = 9)

            if t != tread[-1]:
                dt[k] = float(tread[k+1]) - time[k]
                alphadt = readfield(sol, tread[k+1], 'alpha.a', structured = True, precision = 12)
                Uadt      = readfield(sol, tread[k+1], 'U.a', structured = True, precision = 12)

            else:
                dt[k] = dt[k-1]

            dybdt_l = np.zeros(n1d)
            dybdt_h = np.zeros(n1d)

            dybdx_uphi_l = np.zeros(n1d)
            dybdx_uphi_h = np.zeros(n1d)

            vphi_l = np.zeros(n1d)
            vphi_h = np.zeros(n1d)

            dybdt_phi_l = np.zeros(n1d)
            dybdt_phi_h = np.zeros(n1d)

            Inte_phi   = np.zeros(n1d)
            Inte_phidt = np.zeros(n1d)
            storage    = np.zeros(n1d)

            Inte_uphi = np.zeros(n1d)
            uphiflux  = np.zeros(n1d)

            bed_l   = np.zeros((n1d), dtype = int)
            beddt_l = np.zeros((n1d), dtype = int)
            bed_h   = np.zeros((n1d), dtype = int)
            beddt_h = np.zeros((n1d), dtype = int)

            ybed_l   = np.zeros((n1d))
            ybeddt_l = np.zeros((n1d))
            ybed_h   = np.zeros((n1d))
            ybeddt_h = np.zeros((n1d))

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

            BI   = Budget_Integrand(xbed, Yb, alpha, ybedmethod = 'discrete')
            BIdt = Budget_Integrand(xbed, Yb, alphadt, ybedmethod = 'discrete')

            if xbed[index_xzero]>=0 and xbed[inte]<=xmax - 0.5*ymax:
                xf[k] = np.max(Xb[np.where(alpha>1e-3)])

        ############### find the bed height
                bed_l, bed_h, ybed_l, ybed_h         = BI.bed_height()
                beddt_l, beddt_h, ybeddt_l, ybeddt_h = BIdt.bed_height()

        ########### derivative with space
                dybdx_l = BI.x_deri(ybed_l, xbed)
                dybdx_h = BI.x_deri(ybed_h, xbed)         

        #     print ("bed = \t", ybed_h)

                dybdxdt_l = BIdt.x_deri(ybeddt_l, xbed)
                dybdxdt_h = BIdt.x_deri(ybeddt_h, xbed)

        ########## derivative with time

                dybdt_l = (ybeddt_l - ybed_l) / dt[k]
                dybdt_h = (ybeddt_h - ybed_h) / dt[k]

        #print ("dybdt_h = \t", np.sum(dybdt_h))

                dybeddt_l[:, k] = BI.SurfacevalueMul(alpha[:, :, :], dybdt_l, bed_l)
                dybeddt_h[:, k] = BI.SurfacevalueMul(alpha[:, :, :], dybdt_h, bed_h)

        ######### nonlinear terms

                uphi = alpha*Ua[0, :, :, :]
                vphi = alpha*Ua[1, :, :, :]
        
                uphidt = alphadt*Uadt[0, :, :, :]
                vphidt = alphadt*Uadt[1, :, :, :]

                dybdx_uphi_l = BI.SurfacevalueMul(uphi, dybdx_l, bed_l)
                dybdx_uphi_h = BI.SurfacevalueMul(uphi, dybdx_h, bed_h)

        #print ("dybdx_uphi_l = \t", np.where(dybdx_h<0), "\t dybdx_uphi_h = \t",np.where(dybdx_h>0))

                vphi_l = BI.Surfacevalue(vphi, bed_l)
                vphi_h = BI.Surfacevalue(vphi, bed_h)
        #print ("vphi_l = \t", np.mean(vphi_l))

        ############## Storage

                Inte_phi   = BI.Integrate_flux(alpha, Yb, bed_l, bed_h)
                Inte_phidt = BIdt.Integrate_flux(alphadt, Yb, beddt_l, beddt_h)
                Storage    = (Inte_phidt - Inte_phi) / dt[k]

        #print ("Inte_phi  = \t", np.where(Storage>0), "\t Inte_phi = \t", np.where(Inte_phidt))

        ############# flux
                Inte_phi = BI.Integrate_flux(uphi, Yb, bed_l, bed_h) + BIdt.Integrate_flux(uphidt, Yb, beddt_l, beddt_h)
                uphiflux  = BI.x_deri(Inte_phi, xbed)

        #print ("Inte_phi = \t", np.where())
        ### Right hand side of NS equation
        ############
        ########## Stress value Fluid

        ######## Bouyancy
        ############ Write files

                write_hdf5( direc + '/' + 'Storage_%g.h5'% (time[k]), Storage)
                #write_hdf5('output_2D_%s/Inte_phi_%g.h5'% (time[k]), Inte_phi)
                #write_hdf5('output_2D_%s/Inte_phidt_%g.h5'% (time[k]), Inte_phidt)

        ########### flux
                #write_hdf5('output_2D_%s/dybdx_uphi_l_%g.h5'% (direc,time[k]), dybdx_uphi_l)
                #write_hdf5('output_2D_%s/dybdx_uphi_h_%g.h5'% (direc, time[k]), dybdx_uphi_h)
                #write_hdf5('output_2D_%s/vphi_l_%g.h5'%(direc, time[k]), vphi_l)
                #write_hdf5('output_2D_%s/vphi_h_%g.h5'% (direc, time[k]), vphi_h)
                #write_hdf5('output_2D_%s/uphiflux_%g.h5'% (direc, time[k]), uphiflux)
        ############# bottom part

            else:
                break


        write_hdf5('output_2D_%s/xbed.h5'%(direc), xbed)
        write_hdf5('output_2D_%s/dybeddt_l.h5'%(direc), dybeddt_l)
        write_hdf5('output_2D_%s/dybeddt_h.h5'%(direc), dybeddt_h)


