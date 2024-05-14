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



class Momentum_Budget:
    def __init__(self, path, dim, budget, parameter, value):

        self.path = path
        self.dim = dim
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
        
        #T = T[:5]
        for t in T:
            print ("reading time: %s s"%t)

            k = k + 1

            time[k] = float(t)

#            alpha_a =  readfield(sol, t, 'alpha.a', structured = True, precision = prec )
#            alpha_b =  readfield(sol, t, 'alpha.b', structured = True, precision = prec )
#            Ua      =  fluidfoam.readvector(sol, t, 'U.a', structured = True, precision = prec)
            U       =  fluidfoam.readvector(sol, t, 'U.b', structured = True, precision = prec)
            P       =  fluidfoam.readscalar(sol, t, 'pMech', structured = True, precision = prec)
#            PA      =  fluidfoam.readscalar(sol, t, 'pA', structured = True, precision = prec)
#            tauA    =  fluidfoam.readtensor(sol, t, 'Taua', structured = True, precision = prec)
            tauB    =  fluidfoam.readtensor(sol, t, 'Taub', structured = True, precision = prec)
#            rhoM    =  fluidfoam.readvector(sol, t, 'rho_mixd', structured = True, precision = prec)
#            nonlM   =  fluidfoam.readtensor(sol, t, 'nonl_mixd', structured = True, precision = prec)
#            bouyM   =  fluidfoam.readvector(sol, t, 'bouy_mixd', structured = True, precision = prec)

            if t != tread[-1]:
                dt[k] = float(tread[k+1]) - time[k]
#                alphadt_a = readfield(sol, tread[k+1], 'alpha.a', structured = True, precision = prec)
#                alphadt_b = readfield(sol, tread[k+1], 'alpha.b', structured = True, precision = prec )
#                Uadt      = readfield(sol, tread[k+1], 'U.a', structured = True, precision = prec)
                Udt      = readfield(sol, tread[k+1], 'U.b', structured = True, precision = prec)
                Pdt       = readfield(sol, tread[k+1], 'pMech', structured = True, precision = prec)
#                PAdt      = readfield(sol, tread[k+1], 'pA', structured = True, precision = prec)
#                tauAdt    = readfield(sol, tread[k+1], 'Taua', structured = True, precision = prec)
                tauBdt    = readfield(sol, tread[k+1], 'Taub', structured = True, precision = prec)
#                rhoMdt    = readfield(sol, tread[k+1], 'rho_mixd', structured = True, precision = prec)
#                nonlMdt   = readfield(sol, tread[k+1], 'nonl_mixd', structured = True, precision = prec)
#                bouyMdt   = readfield(sol, tread[k+1], 'bouy_mixd', structured = True, precision = prec)

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
    
            #if int(time[k])==10:

            #   print ("InteA = \t", Inte_A )
            #   print ("IntedtA = \t", Inte_Adt )

            Storage = (Inte_A - Inte_Adt) / dt[k]
    
    ############## flux due to nonlinear term
    
            uxux = np.zeros(n1d)
            uxuxdt = np.zeros(n1d)
            inte_uxux = np.zeros(n1d)
            inte_uxuxdt = np.zeros(n1d)
            flux_uxux = np.zeros(n1d)
    
            uxux   = U[0, :, :, :]* U[0, :, :, :]
            uxuxdt = Udt[0, :, :, :]* Udt[0, :, :, :]
    
            #if int(time[k])==10:
            #   print ("uxux = \t", uxux)
                #  print ("uxuxdt = \t", uxuxdt)
    
            inte_uxux   = BI.Integrate_flux(uxux, Yb, bed_l, bed_h)
            inte_uxuxdt = BIdt.Integrate_flux(uxuxdt, Yb, bed_l, bed_h)
            flux_uxux = (BI.x_deri(inte_uxux, x)  + BIdt.x_deri(inte_uxuxdt, x)) / 2.
    
            #if int(time[k])==10:
            #   print ("flux_uxux = \t", flux_uxux)

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
    
            #if int(time[k])==10:
            #   print ("uxuy_bf = \t", uxuy_bf)
            #  print ("uxuy_hf = \t", uxuy_hf)

    ################# pressure term
    
            xderi_inte_P  = np.zeros(n1d)

            Inte_Press    = BI.Integrate_flux(P[:, :, :], Yb, bed_l, bed_h)
            Intedt_Press  = BIdt.Integrate_flux(Pdt[:, :, :], Yb, beddt_l, beddt_h)

            #if int(time[k])==10:
            #    print ("P = \t", P)
            #    print ("Pdt = \t", Pdt)
            #    print ("Inte P = \t", Inte_Press)
            #    print ("Inte Pdt = \t", Intedt_Press)

    
            Inte_Press_total = (Inte_Press + Intedt_Press) / 2.
            xderi_inte_P = BI.x_deri(Inte_Press_total, x)
    
            #if int(time[k])==10:
                #print ("xderi P = \t",xder_inte_P)

    ############### other technique for pressure

            inte_xderi_P = np.zeros(n1d)

            xder_Press   = BI.x_deri(P[:, :, :], x)
            xder_Pressdt = BIdt.x_deri(Pdt[:, :, :], x)

            inte_xderi_P = (BI.Integrate_flux(xder_Press, Yb, bed_l, bed_h) + BIdt.Integrate_flux(xder_Pressdt, Yb, beddt_l, beddt_h))/2.

    ############ Stress calculations
    
            deri_tau_xx = np.zeros(n1d)
            inte_tau_xx_total = np.zeros(n1d)

    
            inte_tau_xx_total = (BI.Integrate_flux(tauB[0, :, :, :], Yb, bed_l, bed_h) + BIdt.Integrate_flux(tauBdt[0, :, :, :], Yb, bed_l, bed_h)) / 2.
    
            #if int(time[k])==10:
                #print ("inte tau_xx = \t", inte_tau_xx_total)
                #print ("tauB = \t", tauB[0, :, :, :])
            #  print ("tauBdt = \t", tauBdt[0, :, :, :])

            deri_tau_xx = BI.x_deri(inte_tau_xx_total, x)
    
            #if int(time[k])==10:
                #print ("deri tau xx = \t", deri_tau_xx)
    
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
            
            #if int(time[k]==2):
            #    print ("tauxy_b = \t", bed_l)
#                print ("tauxy_h = \t", bed_h)
                #print ("tauxy_b = \t", tauxy_b)

            tauxy_bf = (tauxy_b + tauxydt_b)/2.
            tauxy_hf = (tauxy_h + tauxydt_h)/2.

            #if int(time[k])==10:
     #          print ("tauxy b = \t", tauxy_b)
      #         print ("tauxy h = \t", tauxy_h)
       #        print ("tauxydt  b = \t", tauxydt_b)
        #       print ("tauxydt  h = \t", tauxydt_h)

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
        #write_hdf5("output_2D/y.h5", y)





'''
class Momentum_Budget:
    def __init__(self, path, dire, parameter, value):

        self.path = path
        self.dire = dire
        self.parameter = parameter
        self.value = value

        BF = filter_snapshot(path)

        T       = BF.tprocess
        xbed    = BF.xbed
        Yb      = BF.Yb
        tread   = BF.tread
        Xb = BF.Xb
        ny = BF.ny

        n1d = len(xbed)
        Nt = len(T)
        
        print ("Manohar in Momentum budget class = ##############")

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

        sol = self.path
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

            alpha_a =  readfield(sol, t, 'alpha.a', structured = True, precision = 9 )
            alpha_b =  readfield(sol, t, 'alpha.b', structured = True, precision = 9 )
            Ua      =  fluidfoam.readvector(sol, t, 'U.a', structured = True, precision = 9)
            Ub      =  fluidfoam.readvector(sol, t, 'U.b', structured = True, precision = 9)
            P       =  fluidfoam.readscalar(sol, t, 'pMech', structured = True, precision = 9)
            PA      =  fluidfoam.readscalar(sol, t, 'pA', structured = True, precision = 9)
            tauA    =  fluidfoam.readtensor(sol, t, 'Taua', structured = True, precision = 9)
            tauB    =  fluidfoam.readtensor(sol, t, 'Taub', structured = True, precision = 9)
            rhoM    =  fluidfoam.readvector(sol, t, 'rho_mixd', structured = True, precision = 9)
            nonlM   =  fluidfoam.readtensor(sol, t, 'nonl_mixd', structured = True, precision = 9)
            bouyM   =  fluidfoam.readvector(sol, t, 'bouy_mixd', structured = True, precision = 9)

            if t != tread[-1]:
                dt[k] = float(tread[k+1]) - time[k]
                alphadt_a = readfield(sol, tread[k+1], 'alpha.a', structured = True, precision = 12)
                alphadt_b = readfield(sol, tread[k+1], 'alpha.b', structured = True, precision = 12 )
                Uadt      = readfield(sol, tread[k+1], 'U.a', structured = True, precision = 12)
                Ubdt      = readfield(sol, tread[k+1], 'U.b', structured = True, precision = 12)
                Pdt       = readfield(sol, tread[k+1], 'pMech', structured = True, precision = 12)
                PAdt      = readfield(sol, tread[k+1], 'pA', structured = True, precision = 12)
                tauAdt    = readfield(sol, tread[k+1], 'Taua', structured = True, precision = 12)
                tauBdt    = readfield(sol, tread[k+1], 'Taub', structured = True, precision = 12)
                rhoMdt    = readfield(sol, tread[k+1], 'rho_mixd', structured = True, precision = 12)
                nonlMdt   = readfield(sol, tread[k+1], 'nonl_mixd', structured = True, precision = 12)
                bouyMdt   = readfield(sol, tread[k+1], 'bouy_mixd', structured = True, precision = 12)

            else:
                dt[k] = dt[k-1]

    
            Adybdt_l = np.zeros(n1d)
            Adybdt_h = np.zeros(n1d)
    
            Adybdt_ldt = np.zeros(n1d)
            Adybdt_hdt = np.zeros(n1d)

            dybdt_l = np.zeros(n1d)
            dybdt_h = np.zeros(n1d)

            B_l = np.zeros(n1d)
            B_h = np.zeros(n1d)

            Inte_B = np.zeros(n1d)
            Bflux = np.zeros(n1d)

            A_l   = np.zeros((n1d))
            A_h   = np.zeros((n1d))

            Adt_l   = np.zeros((n1d))
            Adt_h   = np.zeros((n1d))

            B   = np.zeros(n1d)
            Bdt = np.zeros(n1d)

            Acol   = np.zeros(n1d)
            Acoldt = np.zeros(n1d)


            dybdx_B_l   = np.zeros((n1d))
            dybdx_B_h   = np.zeros((n1d))

            C_l   = np.zeros((n1d))
            C_h   = np.zeros((n1d))

            Inte_B = np.zeros(n1d)
            Bflux = np.zeros(n1d)

            Inte_A   = np.zeros(n1d)
            Inte_Adt = np.zeros(n1d)
            Storage  = np.zeros(n1d)

            FPres_l = np.zeros(n1d)
            FPres_h = np.zeros(n1d)
    
            alpha_xderi_FPres = np.zeros((n1d, ny))
            xderi_FPres = np.zeros((n1d, ny))
            xderi_PPres = np.zeros((n1d, ny))

            alphadt_xderi_FPres = np.zeros((n1d, ny))
            xderi_FPresdt = np.zeros((n1d, ny))
            xderi_PPresdt = np.zeros((n1d, ny))

            #FPres_h = np.zeros(n1d)
    
            FInte_Pres = np.zeros(n1d)
            FPresflux = np.zeros(n1d)

            PPres_l = np.zeros(n1d)
            PPres_h = np.zeros(n1d)
            PInte_Pres = np.zeros(n1d)
            PPresflux = np.zeros(n1d)

            dybdx_TauP_l = np.zeros(n1d)
            dybdx_TauP_h = np.zeros(n1d)
            Inte_TauP = np.zeros(n1d)
            TauPflux = np.zeros(n1d)

            dybdx_TauF_l = np.zeros(n1d)
            dybdx_TauF_h = np.zeros(n1d)
            Inte_TauF = np.zeros(n1d)
            TauFflux = np.zeros(n1d)

            TauxyF_l = np.zeros(n1d)
            TauxyF_h = np.zeros(n1d)

            TauxyP_l = np.zeros(n1d)
            TauxyP_h = np.zeros(n1d)
    

            Bouy = np.zeros(n1d)

            bed_l   = np.zeros((n1d), dtype = int)
            beddt_l = np.zeros((n1d), dtype = int)
            bed_h   = np.zeros((n1d), dtype = int)
            beddt_h = np.zeros((n1d), dtype = int)
    
            ybed_l   = np.zeros((n1d))
            ybeddt_l = np.zeros((n1d))
            ybed_h   = np.zeros((n1d))
            ybeddt_h = np.zeros((n1d))

            dybdx_l   = np.zeros((n1d))
            dybdxdt_l = np.zeros((n1d))
            dybdx_h   = np.zeros((n1d))
            dybdxdt_h = np.zeros((n1d))
    
            dybdxa   = np.zeros((n1d))
            dybdxadt = np.zeros((n1d))
    

            xf[k] = np.max(Xb[np.where(alpha_a>1e-3)])
            range_inte = np.where(alpha_a[:, 0, 0]>1e-3)[0]
            inte = range_inte[-1]

            range_inte_dt = np.where(alpha_a[:, 0, 0]>1e-3)[0]
            inte_dt = range_inte_dt[-1]

            xzero = np.where(xbed>0)[0]
            index_xzero =  xzero[0]

            print ("print xzero = \t", index_xzero)
            print ("print inte = \t", inte)

            f.write('%d \t %g \t %g \n'%(time[k], index_xzero, inte))

            xf[k] = np.max(Xb[np.where(alpha_a>1e-3)])

            BI   = Budget_Integrand(xbed, Yb, alpha_a, ybedmethod = 'discrete')
            BIdt = Budget_Integrand(xbed, Yb, alphadt_a, ybedmethod = 'discrete')

            if xbed[index_xzero]>=0 and xbed[inte]<=xmax - 0.5*ymax:

    ############### find the bed height
                bed_l, bed_h, ybed_l, ybed_h = BI.bed_height()
                beddt_l, beddt_h, ybeddt_l, ybeddt_h = BIdt.bed_height()
    
    ########### derivative with space
                dybdx_l = BI.x_deri(ybed_l, xbed)         
                dybdx_h = BI.x_deri(ybed_h, xbed)         

                dybdxdt_l = BIdt.x_deri(ybeddt_l, xbed)         
                dybdxdt_h = BIdt.x_deri(ybeddt_h, xbed)  

    #    print ("shape of rhoM = \t", nonlM.shape)
        
        ########## derivative with time

                dybdt_l = (ybeddt_l - ybed_l) / dt[k]
                dybdt_h = (ybeddt_h - ybed_h) / dt[k]


                dybeddt_l[:, k] = BI.SurfacevalueMul(rhoM[0, :, :, :], dybdt_l, bed_l) 
                dybeddt_h[:, k] = BI.SurfacevalueMul(rhoM[0, :, :, :], dybdt_h, bed_h) 
    
        ######### nonlinear terms

                B_l = (BI.SurfacevalueMul(nonlM[0, :, :, :], dybdx_l, bed_l) + BIdt.SurfacevalueMul(nonlMdt[0, :, :, :], dybdxdt_l, beddt_l))/2.
                B_h = (BI.SurfacevalueMul(nonlM[0, :, :, :], dybdx_h, bed_h) + BIdt.SurfacevalueMul(nonlMdt[0, :, :, :], dybdxdt_h, beddt_h))/2.
                #B_h = SurfacevalueMul(nonlM[0, :, :, :], dybdx_h, bed_h)
        
                C_l = (BI.SurfacevalueMul(nonlM[1, :, :, :], dybdx_l, bed_l) + BIdt.SurfacevalueMul(nonlMdt[1, :, :, :], dybdxdt_l, beddt_l))/2.
                C_h = (BI.SurfacevalueMul(nonlM[1, :, :, :], dybdx_h, bed_h) + BIdt.SurfacevalueMul(nonlMdt[1, :, :, :], dybdxdt_h, beddt_h))/2.
    
                #C_l = Surfacevalue(nonlM[1, :, :, :], bed_l)
                #C_h = Surfacevalue(nonlM[1, :, :, :], bed_h)

        ############## Storage
    
                Acol   = rhoM[0, :, :, :]    #rho_a* alpha_a*Ua[0, :, :, :]
                Acoldt = rhoMdt[0, :, :, :]    #rho_a* alphadt_a*Uadt[0, :, :, :]
    
                dybeddt_l[:, k] = BI.SurfacevalueMul(Acol, dybdt_l, bed_l) 
                dybeddt_h[:, k] = BIdt.SurfacevalueMul(Acoldt, dybdt_h, bed_h) 

                Inte_A   = BI.Integrate_flux(Acol, Yb, bed_l, bed_h)
                Inte_Adt = BIdt.Integrate_flux(Acoldt, Yb, beddt_l, beddt_h)
                Storage  = (Inte_Adt - Inte_A) / dt[k] 

        ############# flux
    
                B = nonlM[0, :, :, :]     #rho_a* alpha_a* Ua[0, :, :, :]* Ua[0, :, :, :]
                Bdt = nonlMdt[0, :, :, :] #rho_a* alphadt_a* Uadt[0, :, :, :]* Uadt[0, :, :, :]

                Inte_B   = BI.Integrate_flux(B, Yb, bed_l, bed_h)
                #print ("Inte_B = \t", Inte_B)
                Inte_Bdt = BIdt.Integrate_flux(Bdt, Yb, beddt_l, beddt_h) 
                Bflux  = (BI.x_deri(Inte_B, xbed) + BIdt.x_deri(Inte_Bdt, xbed))/2.

        ### Right hand side of NS equation
        ## Fluid pressure
                FPres_l    = (BI.SurfacevalueMul(P, dybdx_l, bed_l) + BIdt.SurfacevalueMul(Pdt, dybdxdt_l, beddt_l))/2.
                FPres_h    = (BI.SurfacevalueMul(P, dybdx_h, bed_h) + BIdt.SurfacevalueMul(Pdt, dybdxdt_h, beddt_h))/2.
    
                FInte_Pres = (BI.Integrate_flux(P, Yb, bed_l, bed_h) + BIdt.Integrate_flux(Pdt, Yb, beddt_l, beddt_h)) / 2.
                FPressflux = BI.x_deri(FInte_Pres, xbed)

        ## Particle pressure
                PPres_l     = (BI.SurfacevalueMul(PA, dybdx_l, bed_l) + BIdt.SurfacevalueMul(PAdt, dybdxdt_l, beddt_h))/2.
                PPres_h     = (BI.SurfacevalueMul(PA, dybdx_h, bed_h) + BIdt.SurfacevalueMul(PAdt, dybdxdt_h, beddt_h))/2.
                PInte_Pres  = (BI.Integrate_flux(PA, Yb, bed_l, bed_h) + BIdt.Integrate_flux(PAdt, Yb, beddt_l, beddt_h))/2.
                PPressflux =  BI.x_deri(PInte_Pres, xbed)
        ########## Stress value particle

                dybdx_TauP_l = (BI.SurfacevalueMul(tauA[0, :, :, :], dybdx_l, bed_l) + BIdt.SurfacevalueMul(tauAdt[0, :, :, :], dybdxdt_l, beddt_l))/2.
                dybdx_TauP_h = (BI.SurfacevalueMul(tauA[0, :, :, :], dybdx_h, bed_h) + BIdt.SurfacevalueMul(tauAdt[0, :, :, :], dybdxdt_h, beddt_h))/2.
                Inte_TauP    = (BI.Integrate_flux(tauA[0, :, :, :], Yb, bed_l, bed_h) + BIdt.Integrate_flux(tauAdt[0, :, :, :], Yb, beddt_l, beddt_h))/2.
                TauPflux     = BI.x_deri(Inte_TauP, xbed)
                TauxyP_l     = (BI.Surfacevalue(tauA[1, :, :, : ], bed_l) + BIdt.Surfacevalue(tauAdt[1, :, :, :], beddt_l))/2.
                TauxyP_h     = (BI.Surfacevalue(tauA[1, :, :, : ], bed_h) + BIdt.Surfacevalue(tauAdt[1, :, :, :], beddt_h))/2.
                #TauxyP_h     = Surfacevalue(tauA[1, :, :, : ], bed_h)
        
    
        ############
        ########## Stress value Fluid
                dybdx_TauF_l = (BI.SurfacevalueMul(tauB[0, :, :, :], dybdx_l, bed_l) + BIdt.SurfacevalueMul(tauBdt[0, :, :, :], dybdxdt_l, beddt_l))/2.
                dybdx_TauF_h = (BI.SurfacevalueMul(tauB[0, :, :, :], dybdx_h, bed_h) + BIdt.SurfacevalueMul(tauBdt[0, :, :, :], dybdxdt_h, beddt_h))/2.
                Inte_TauF    = (BI.Integrate_flux(tauB[0, :, :, :], Yb, bed_l, bed_h) + BIdt.Integrate_flux(tauBdt[0, :, :, :], Yb, beddt_l, beddt_h))/2.
                TauFflux     = BI.x_deri(Inte_TauF, xbed)
                TauxyF_l     = (BI.Surfacevalue(tauB[1, :, :, : ], bed_l) + BIdt.Surfacevalue(tauBdt[1, :, :, :], beddt_l))/2.
                TauxyF_h     = (BI.Surfacevalue(tauB[1, :, :, : ], bed_h) + BIdt.Surfacevalue(tauBdt[1, :, :, :], beddt_h))/2.
                #TauxyP_h     = Surfacevalue(tauA[1, :, :, : ], bed_h)
       
        ####################
                gx = 0
                gy = 0
                #Bouy = (Integrate_flux( rho_a* alpha_a*gx , Yb, bed_l, bed_h) + Integrate_flux(rho_a*alphadt_a*gx, Yb, beddt_l, beddt_h))/2.
                Bouy = (BI.Integrate_flux( bouyM[0, :, :, :] , Yb, bed_l, bed_h) + BIdt.Integrate_flux(bouyMdt[0, :, :, :], Yb, beddt_l, beddt_h))/2.

        ############ Write files

                write_hdf5( direc + '/' + 'Storage_%g.h5'%(time[k]), Storage)

        ########### flux
#                write_hdf5('output_2D/Bflux_%g.h5'%(time[k]), Bflux)
#                write_hdf5('output_2D/FPresflux_%g.h5'%(time[k]), FPressflux)
#                write_hdf5('output_2D/PPresflux_%g.h5'%(time[k]), PPressflux)
#                write_hdf5('output_2D/TauFflux_%g.h5'%(time[k]), TauFflux)
#                write_hdf5('output_2D/TauPflux_%g.h5'%(time[k]), TauPflux)
       
        ############# bottom part
        
#                write_hdf5('output_2D/B_l_%g.h5'%(time[k]), B_l)
#                write_hdf5('output_2D/C_l_%g.h5'%(time[k]), C_l)
#                write_hdf5('output_2D/TauFxy_l_%g.h5'%(time[k]), TauxyF_l)
#                write_hdf5('output_2D/TauPxy_l_%g.h5'%(time[k]), TauxyP_l)
#                write_hdf5('output_2D/dybdx_TauFxy_l_%g.h5'%(time[k]), dybdx_TauF_l)
#                write_hdf5('output_2D/dybdx_TauPxy_l_%g.h5'%(time[k]), dybdx_TauP_l)
#                write_hdf5('output_2D/fluidP_l_%g.h5'%(time[k]), FPres_l)
#                write_hdf5('output_2D/particleP_l_%g.h5'%(time[k]), PPres_l)
        
        ############# top part
#                write_hdf5('output_2D/C_h_%g.h5'%(time[k]), C_h)
#                write_hdf5('output_2D/TauFxy_h_%g.h5'%(time[k]), TauxyF_h)
#                write_hdf5('output_2D/TauPxy_h_%g.h5'%(time[k]), TauxyP_h)
#                write_hdf5('output_2D/dybdx_TauFxy_h_%g.h5'%(time[k]), dybdx_TauF_h)
#                write_hdf5('output_2D/dybdx_TauPxy_h_%g.h5'%(time[k]), dybdx_TauP_h)
#                write_hdf5('output_2D/fluidP_h_%g.h5'%(time[k]), FPres_l)
#                write_hdf5('output_2D/particleP_h_%g.h5'%(time[k]), PPres_h)
        
        ################ Bouyancy

#                write_hdf5('output_2D/bouy_%g.h5'%(time[k]), Bouy)

            else:
                break


#        write_hdf5('output_2D/xbed.h5', xbed)
#        write_hdf5('output_2D/dybeddt_l.h5', dybeddt_l)
#        write_hdf5('output_2D/dybeddt_h.h5', dybeddt_h)
'''




