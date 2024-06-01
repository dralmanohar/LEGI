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

from Operations.spatial_operation.two.Budget_Ingradient import Budget_Integrand
from Operations.read_write.diagnostic_function import write_hdf5
from Operations.extract_time.Filter_Snapshot import Filter_Snapshot



class Momentum_Budget_Poiseuille:
    def __init__(self, path, dim, budget, parameter, value):

        self.path = path
        self.dim = dimz
        self.budget = budget
        self.parameter = parameter
        self.value = value
        
        sol = self.path

        BF = Filter_Snapshot(sol, self.dim)

        T       = BF.tprocess
        x       = BF.xbed
        Yb      = BF.Yb
        tread   = BF.tread
        Xb = BF.Xb
        ny = BF.ny

        n1d = len(x)
        Nt = len(T)

        print ("n1d = \t", n1d)
        print ("Nt = \t", Nt)

        time        = np.zeros(Nt)
        dt          = np.zeros(Nt)
        dybeddt_l   = np.zeros((n1d, Nt))
        dybeddt_h   = np.zeros((n1d, Nt))
        balance     = np.zeros((n1d, Nt))
        xf          = np.zeros(Nt)
        
        xbed = x
        #xmax = np.amax(xbed)
        #ymax = np.amax(Yb[0, :, 0])

        k = -1

        prec = 9
        #f = open('index_zero_xf.d', 'w')

        direc = './output_%s_%s/%s/%s'%(self.dim, self.budget, self.parameter, self.value)

        try:
            os.makedirs(direc, exist_ok = True)
            print("Directory '%s' created successfully" % direc)
        except OSError as error:
            print("Directory '%s' can not be created" % direc)
        
        T = T[:10]
        for t in T:
            print ("reading time: %s s"%t)

            k = k + 1

            time[k] = float(t)

            alpha_a =  readfield(sol, t, 'alpha.a', structured = True, precision = prec )
            alpha_b =  readfield(sol, t, 'alpha.b', structured = True, precision = prec )
            Ua      =  fluidfoam.readvector(sol, t, 'U.a', structured = True, precision = prec)
            U       =  fluidfoam.readvector(sol, t, 'U.b', structured = True, precision = prec)
            P       =  fluidfoam.readscalar(sol, t, 'pMech', structured = True, precision = prec)
            PA      =  fluidfoam.readscalar(sol, t, 'pA', structured = True, precision = prec)
            tauA    =  fluidfoam.readtensor(sol, t, 'Taua', structured = True, precision = prec)
            tauB    =  fluidfoam.readtensor(sol, t, 'Taub', structured = True, precision = prec)
            rhoM    =  fluidfoam.readvector(sol, t, 'rho_mixd', structured = True, precision = prec)
            nonlM   =  fluidfoam.readtensor(sol, t, 'nonl_mixd', structured = True, precision = prec)
            bouyM   =  fluidfoam.readvector(sol, t, 'bouy_mixd', structured = True, precision = prec)

            if t != tread[-1]:
                dt[k] = float(tread[k+1]) - time[k]
                alphadt_a = readfield(sol, tread[k+1], 'alpha.a', structured = True, precision = prec)
                alphadt_b = readfield(sol, tread[k+1], 'alpha.b', structured = True, precision = prec )
                Uadt      = readfield(sol, tread[k+1], 'U.a', structured = True, precision = prec)
                Udt      = readfield(sol, tread[k+1], 'U.b', structured = True, precision = prec)
                Pdt       = readfield(sol, tread[k+1], 'pMech', structured = True, precision = prec)
                PAdt      = readfield(sol, tread[k+1], 'pA', structured = True, precision = prec)
                tauAdt    = readfield(sol, tread[k+1], 'Taua', structured = True, precision = prec)
                tauBdt    = readfield(sol, tread[k+1], 'Taub', structured = True, precision = prec)
                rhoMdt    = readfield(sol, tread[k+1], 'rho_mixd', structured = True, precision = prec)
                nonlMdt   = readfield(sol, tread[k+1], 'nonl_mixd', structured = True, precision = prec)
                bouyMdt   = readfield(sol, tread[k+1], 'bouy_mixd', structured = True, precision = prec)

            else:
                dt[k] = dt[k-1]

            ########## define the instance 

            BI   = Budget_Integrand(x, Yb)#Budget_Integrand(x, Yb)
            BIdt = Budget_Integrand(x, Yb)#Budget_Integrand(x, Yb)

            ################## calculation of bed height

            ####### calculation of bed height

            bed_l   = np.zeros((n1d), dtype = int)
            beddt_l = np.zeros((n1d), dtype = int)
            bed_h   = np.zeros((n1d), dtype = int)
            beddt_h = np.zeros((n1d), dtype = int)

            ybed_l   = np.zeros((n1d))
            ybeddt_l = np.zeros((n1d))
            ybed_h   = np.zeros((n1d))
            ybeddt_h = np.zeros((n1d))

            bed_l, bed_h, ybed_l, ybed_h = BI.bed_height()
            beddt_l, beddt_h, ybeddt_l, ybeddt_h = BIdt.bed_height()

            ########### derivative with space

            dybdx_l = np.zeros(n1d)
            dybdx_h = np.zeros(n1d)

            dybdxdt_l = np.zeros(n1d)
            dybdxdt_h = np.zeros(n1d)

            dybdx_l = BI.x_deri(ybed_l, xbed)
            dybdx_h = BI.x_deri(ybed_h, xbed)

            dybdxdt_l = BIdt.x_deri(ybeddt_l, xbed)
            dybdxdt_h = BIdt.x_deri(ybeddt_h, xbed)
            

            ########## derivative with time

            dybdt_l = np.zeros(n1d)
            dybdt_h = np.zeros(n1d)
            
            dybddt_l = np.zeros(n1d)
            dybddt_h = np.zeros(n1d)


            dybdt_l = (ybeddt_l - ybed_l) / dt[k]
            dybdt_h = (ybeddt_h - ybed_h) / dt[k]

            dybeddt_l[:, k] = BI.SurfacevalueMul(rhoM[0, :, :, :], dybdt_l, bed_l)
            dybeddt_h[:, k] = BI.SurfacevalueMul(rhoM[0, :, :, :], dybdt_h, bed_h) 



            ##### Storage terms 

            Acol   = np.zeros(n1d)
            Acoldt = np.zeros(n1d)

            Inte_A   = np.zeros(n1d)
            Inte_Adt = np.zeros(n1d)
            Storage  = np.zeros(n1d)

            Acol = U[0, :, :, :]
            Acoldt = Udt[0, :, :, :]
            
        
            Inte_A   = BI.Integrate_flux(Acol, Yb, bed_l, bed_h)
            Inte_Adt = BIdt.Integrate_flux(Acoldt, Yb, beddt_l, beddt_h)
            Storage = (Inte_Adt - Inte_A) / dt[k]

            ############## flux due to nonlinear term
        
            uxux = np.zeros(n1d)
            uxuxdt = np.zeros(n1d)
            inte_uxux = np.zeros(n1d)
            inte_uxuxdt = np.zeros(n1d)
            flux_uxux = np.zeros(n1d)
    
            uxux   = U[0, :, :, :]* U[0, :, :, :]
            uxuxdt = Udt[0, :, :, :]* Udt[0, :, :, :]
            
    
            inte_uxux   = BI.Integrate_flux(uxux, Yb, bed_l, bed_h)
            inte_uxuxdt = BIdt.Integrate_flux(uxuxdt, Yb, bed_l, bed_h)
            flux_uxux = (BI.x_deri(inte_uxux, x)  + BIdt.x_deri(inte_uxuxdt, x)) / 2.
    
            ################################## Stress value particle


            ############ Variables at the top and bottom of control volume
    
            uxuy   = U[0, :, :, :]*U[1, :, :, :]
            uxuydt = Udt[0, :, :, :]*Udt[1, :, :, :]

            uxuy_b = np.zeros(n1d)
            uxuy_h = np.zeros(n1d)
    
            uxuydt_b = np.zeros(n1d)
            uxuydt_h = np.zeros(n1d)

            uxuy_bf = np.zeros(n1d)
            uxuy_hf = np.zeros(n1d)
    

            uxuy_b = BI.Surfacevalue(uxuy, bed_l)
            uxuy_h = BI.Surfacevalue(uxuy, bed_h)

            uxuydt_b = BIdt.Surfacevalue(uxuydt, beddt_l)
            uxuydt_h = BIdt.Surfacevalue(uxuydt, beddt_h)

            uxuy_bf = (uxuy_b + uxuydt_b)/2.
            uxuy_hf = (uxuy_h + uxuydt_h)/2.
     
   
            ############## pressure derivative

            xderi_inte_P  = np.zeros(n1d)

            Inte_Press    = BI.Integrate_flux(P[:, :, :], Yb, bed_l, bed_h)
            Intedt_Press  = BIdt.Integrate_flux(Pdt[:, :, :], Yb, beddt_l, beddt_h)
    
            Inte_Press_total = (Inte_Press + Intedt_Press) / 2.
            xderi_inte_P = BI.x_deri(Inte_Press_total, x)

            ############### other technique for pressure

            inte_xderi_P = np.zeros(n1d)
            xder_Press   = BI.x_deri(P[:, :, :], x)
            xder_Pressdt = BIdt.x_deri(Pdt[:, :, :], x)
            inte_xderi_P = (BI.Integrate_flux(xder_Press, Yb, bed_l, bed_h) + BIdt.Integrate_flux(xder_Pressdt, Yb, beddt_l, beddt_h))/2.

            ############ Stress calculations

            deri_tau_xx = np.zeros(n1d)
            inte_tau_xx_total = np.zeros(n1d)

            inte_tau_xx_total = (BI.Integrate_flux(tauB[0, :, :, :], Yb, bed_l, bed_h) + BIdt.Integrate_flux(tauBdt[0, :, :, :], Yb, bed_l, bed_h)) / 2.

            deri_tau_xx = BI.x_deri(inte_tau_xx_total, x)

            ################ Stress at the top and bottom

            tauxy_b = np.zeros(n1d)
            tauxy_h = np.zeros(n1d)
            tauxydt_b = np.zeros(n1d)
            tauxydt_h = np.zeros(n1d)

            tauxy_bf = np.zeros(n1d)
            tauxy_hf = np.zeros(n1d)

            tauxy_b = BI.Surfacevalue(tauB[1, :, :, :], bed_l)
            tauxy_h = BI.Surfacevalue(tauB[1, :, :, :], bed_h)

            tauxydt_b = BIdt.Surfacevalue(tauBdt[1, :, :, :], beddt_l)
            tauxydt_h = BIdt.Surfacevalue(tauBdt[1, :, :, :], beddt_h)

            tauxy_bf = (tauxy_b + tauxydt_b)/2.
            tauxy_hf = (tauxy_h + tauxydt_h)/2.


        ############### writing the variables in the hdf5 file
            write_hdf5(direc + '/' + "Storage_%g.h5"%(time[k]), Storage)
            write_hdf5(direc + '/' + "flux_xx_%g.h5"%(time[k]), flux_uxux)
            write_hdf5(direc + '/' + "uxuy_b_%g.h5"%(time[k]), uxuy_bf)
            write_hdf5(direc + '/' + "uxuy_h_%g.h5"%(time[k]), uxuy_hf)
            write_hdf5(direc + '/' + "xderi_inte_P_%g.h5"%(time[k]), xderi_inte_P)
            write_hdf5(direc + '/' + "inte_xderi_P_%g.h5"%(time[k]), inte_xderi_P)
            write_hdf5(direc + '/' + "xderi_tau_xx_%g.h5"%(time[k]), deri_tau_xx)
            write_hdf5(direc + '/' + "tau_xy_b_%g.h5"%(time[k]), tauxy_bf)
            write_hdf5(direc + '/' + "tau_xy_h_%g.h5"%(time[k]), tauxy_hf)

        write_hdf5(direc + '/' + "xbed.h5", x)





