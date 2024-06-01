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

from Operations.spatial_operation.three.Budget_Ingradient import Budget_Integrand
from Operations.read_write.diagnostic_function import write_hdf5
from Operations.extract_time.Filter_Snapshot import Filter_Snapshot



class Momentum_Budget_Turbidity:
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
        Zb = BF.zbed
        ny = BF.ny

        xbed = x

        n1d, n2d = xbed.shape #len(x)
        Nt = len(T)


        time        = np.zeros(Nt)
        dt          = np.zeros(Nt)
        dybeddt_l   = np.zeros((n1d, n2d, Nt))
        dybeddt_h   = np.zeros((n1d, n2d, Nt))
        balance     = np.zeros((n1d, n2d, Nt))
        xf          = np.zeros(Nt)
        
        xbed = x
        xmax = np.amax(np.mean(xbed, axis = 1))
        ymax = np.amax(Yb[0, :, 0])
        k = -1

        prec = 9
        f = open('index_zero_xf.d', 'w')

        direc = './output_%s_%s/%s/%s'%(self.dim, self.budget, self.parameter, self.value)

        try:
            os.makedirs(direc, exist_ok = True)
            print("Directory '%s' created successfully" % direc)
        except OSError as error:
            print("Directory '%s' can not be created" % direc)
        
        T = T[10:14]

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

            ########################

            xf[k] = np.mean(np.max(Xb[np.where(alpha_a>1e-3)]))


            range_inte = np.where(alpha_a[:, 0, 0]>1e-3)[0]
            inte = range_inte[-1]

            range_inte_dt = np.where(alpha_a[:, 0, 0]>1e-3)[0]
            inte_dt = range_inte_dt[-1]

            xzero = np.where(xbed>0)[0]
            index_xzero =  xzero[0]


            f.write('%d \t %g \n'%(time[k], xf[k]))

            BI   = Budget_Integrand(xbed, Yb, alpha=alpha_a)
            BIdt = Budget_Integrand(xbed, Yb, alpha=alphadt_a)
                
            if xf[k]<=xmax - 0.5*ymax:
                #xf[k] = np.max(Xb[np.where(alpha>1e-3)])
                print ("Manohr in the loop")

                ########## define the instance 

                BI   = Budget_Integrand(x, Yb)#Budget_Integrand(x, Yb)
                BIdt = Budget_Integrand(x, Yb)#Budget_Integrand(x, Yb)

                ################## calculation of bed height
                ####### calculation of bed height

                bed_l   = np.zeros((n1d, n2d), dtype = int)
                beddt_l = np.zeros((n1d, n2d), dtype = int)
                bed_h   = np.zeros((n1d, n2d), dtype = int)
                beddt_h = np.zeros((n1d, n2d), dtype = int)

                ybed_l   = np.zeros((n1d, n2d))
                ybeddt_l = np.zeros((n1d, n2d))
                ybed_h   = np.zeros((n1d, n2d))

                ybeddt_h = np.zeros((n1d, n2d))

                bed_l, bed_h, ybed_l, ybed_h = BI.bed_height()
                beddt_l, beddt_h, ybeddt_l, ybeddt_h = BIdt.bed_height()


                ########### derivative with space

                dybdx_l = np.zeros((n1d, n2d))
                dybdx_h = np.zeros((n1d, n2d))

                dybdxdt_l = np.zeros((n1d, n2d))
                dybdxdt_h = np.zeros((n1d, n2d))

                dybdx_l = BI.x_deri(ybed_l, xbed)
                dybdx_h = BI.x_deri(ybed_h, xbed)

                dybdxdt_l = BIdt.x_deri(ybeddt_l, xbed)
                dybdxdt_h = BIdt.x_deri(ybeddt_h, xbed)


                ########## derivative with time

                dybdt_l = np.zeros((n1d, n2d))
                dybdt_h = np.zeros((n1d, n2d))

                dybddt_l = np.zeros((n1d, n2d))
                dybddt_h = np.zeros((n1d, n2d))

                dybdt_l = (ybeddt_l - ybed_l) / dt[k]
                dybdt_h = (ybeddt_h - ybed_h) / dt[k]

                dybeddt_l[:, :, k] = BI.SurfacevalueMul(rhoM[0, :, :, :], dybdt_l, bed_l)
                dybeddt_h[:, :, k] = BI.SurfacevalueMul(rhoM[0, :, :, :], dybdt_h, bed_h) 

                ######### nonlinear terms

                B_l = np.zeros((n1d, n2d))
                B_h = np.zeros((n1d, n2d))
                C_l = np.zeros((n1d, n2d))
                C_h = np.zeros((n1d, n2d))

                B_l = (BI.SurfacevalueMul(nonlM[0, :, :, :], dybdx_l, bed_l) + BIdt.SurfacevalueMul(nonlMdt[0, :, :, :], dybdxdt_l, beddt_l))/2.
                B_h = (BI.SurfacevalueMul(nonlM[0, :, :, :], dybdx_h, bed_h) + BIdt.SurfacevalueMul(nonlMdt[0, :, :, :], dybdxdt_h, beddt_h))/2.

                C_l = (BI.SurfacevalueMul(nonlM[1, :, :, :], dybdx_l, bed_l) + BIdt.SurfacevalueMul(nonlMdt[1, :, :, :], dybdxdt_l, beddt_l))/2.
                C_h = (BI.SurfacevalueMul(nonlM[1, :, :, :], dybdx_h, bed_h) + BIdt.SurfacevalueMul(nonlMdt[1, :, :, :], dybdxdt_h, beddt_h))/2.


                ##### Storage terms 

                Acol   = np.zeros((n1d, n2d))
                Acoldt = np.zeros((n1d, n2d))

                Inte_A   = np.zeros((n1d, n2d))
                Inte_Adt = np.zeros((n1d, n2d))
                Storage  = np.zeros((n1d, n2d))


                Acol = rhoM[0, :, :, :]
                Acoldt = rhoMdt[0, :, :, :]


                Inte_A   = BI.Integrate_flux(Acol, Yb, bed_l, bed_h)
                Inte_Adt = BIdt.Integrate_flux(Acoldt, Yb, beddt_l, beddt_h)
                
                Storage = (Inte_Adt - Inte_A) / dt[k]


                ############## flux due to nonlinear term

                Inte_B = np.zeros((n1d, n2d))
                Inte_Bdt = np.zeros((n1d, n2d))
                Bflux = np.zeros((n1d, n2d))

                B = nonlM[0, :, :, :] 
                Bdt = nonlMdt[0, :, :, :] 


                Inte_B   = BI.Integrate_flux(B, Yb, bed_l, bed_h)
                Inte_Bdt = BIdt.Integrate_flux(Bdt, Yb, beddt_l, beddt_h)
                Bflux  = (BI.x_deri(Inte_B, xbed) + BIdt.x_deri(Inte_Bdt, xbed))/2.

                print ("Bflux = \t", Bflux)


                ###################################### Pressure for fluid

                FPres_l = np.zeros((n1d, n2d))
                FPres_h = np.zeros((n1d, n2d))

                FInte_Pres = np.zeros((n1d, n2d))
                FPressflux = np.zeros((n1d, n2d))

                FPres_l    = (BI.SurfacevalueMul(P, dybdx_l, bed_l) + BIdt.SurfacevalueMul(Pdt, dybdxdt_l, beddt_l))/2.
                FPres_h    = (BI.SurfacevalueMul(P, dybdx_h, bed_h) + BIdt.SurfacevalueMul(Pdt, dybdxdt_h, beddt_h))/2.

                FInte_Pres = (BI.Integrate_flux(P, Yb, bed_l, bed_h) + BIdt.Integrate_flux(Pdt, Yb, beddt_l, beddt_h)) / 2.
                FPressflux = BI.x_deri(FInte_Pres, xbed)


                ################################## Stress value particle

                dybdx_Taup_l = np.zeros((n1d, n2d))
                dybdx_Taup_h = np.zeros((n1d, n2d))
                Inte_TauP  = np.zeros((n1d, n2d))
                TauPflux  = np.zeros((n1d, n2d))

                TauxyP_l  = np.zeros((n1d, n2d))
                TauxyP_h  = np.zeros((n1d, n2d))

                dybdx_TauP_l = (BI.SurfacevalueMul(tauA[0, :, :, :], dybdx_l, bed_l) + BIdt.SurfacevalueMul(tauAdt[0, :, :, :], dybdxdt_l, beddt_l))/2.
                dybdx_TauP_h = (BI.SurfacevalueMul(tauA[0, :, :, :], dybdx_h, bed_h) + BIdt.SurfacevalueMul(tauAdt[0, :, :, :], dybdxdt_h, beddt_h))/2.
                Inte_TauP    = (BI.Integrate_flux(tauA[0, :, :, :], Yb, bed_l, bed_h)+ BIdt.Integrate_flux(tauAdt[0, :, :, :], Yb, beddt_l, beddt_h))/2.
                TauPflux     = BI.x_deri(Inte_TauP, xbed)
                TauxyP_l     = (BI.Surfacevalue(tauA[1, :, :, : ], bed_l) + BIdt.Surfacevalue(tauAdt[1, :, :, :], beddt_l))/2.
                TauxyP_h     = (BI.Surfacevalue(tauA[1, :, :, : ], bed_h) + BIdt.Surfacevalue(tauAdt[1, :, :, :], beddt_h))/2.

                ###################### Stress value Fluid

                dybdx_TauF_l = np.zeros((n1d, n2d))
                dybdx_TauF_h = np.zeros((n1d, n2d))
                Inte_TauF = np.zeros((n1d, n2d))
                TauFflux = np.zeros((n1d, n2d))
                TauxyF_l = np.zeros((n1d, n2d))
                TauxyF_h = np.zeros((n1d, n2d))

                dybdx_TauF_l = (BI.SurfacevalueMul(tauB[0, :, :, :], dybdx_l, bed_l) + BIdt.SurfacevalueMul(tauBdt[0, :, :, :], dybdxdt_l, beddt_l))/2.
                dybdx_TauF_h = (BI.SurfacevalueMul(tauB[0, :, :, :], dybdx_h, bed_h) + BIdt.SurfacevalueMul(tauBdt[0, :, :, :], dybdxdt_h, beddt_h))/2.
                Inte_TauF    = (BI.Integrate_flux(tauB[0, :, :, :], Yb, bed_l, bed_h) + BIdt.Integrate_flux(tauBdt[0, :, :, :], Yb, beddt_l, beddt_h))/2.
                TauFflux     = BI.x_deri(Inte_TauF, xbed)

                TauxyF_l     = (BI.Surfacevalue(tauB[1, :, :, : ], bed_l) + BIdt.Surfacevalue(tauBdt[1, :, :, :], beddt_l))/2.
                TauxyF_h     = (BI.Surfacevalue(tauB[1, :, :, : ], bed_h) + BIdt.Surfacevalue(tauBdt[1, :, :, :], beddt_h))/2.
                
                ################### particle pressure ## Particle pressure

                PPres_l = np.zeros((n1d, n2d))
                PPres_h = np.zeros((n1d, n2d))
                PInte_Pres = np.zeros((n1d, n2d))
                PPressflux = np.zeros((n1d, n2d))

                PPres_l     = (BI.SurfacevalueMul(PA, dybdx_l, bed_l) + BIdt.SurfacevalueMul(PAdt, dybdxdt_l, beddt_h))/2.
                PPres_h     = (BI.SurfacevalueMul(PA, dybdx_h, bed_h) + BIdt.SurfacevalueMul(PAdt, dybdxdt_h, beddt_h))/2.
                PInte_Pres  = (BI.Integrate_flux(PA, Yb, bed_l, bed_h) + BIdt.Integrate_flux(PAdt, Yb, beddt_l, beddt_h))/2.
                PPressflux =  BI.x_deri(PInte_Pres, xbed)


                ############ Variables at the top and bottom of control volume
                ############### writing the variables in the hdf5 file

                write_hdf5(direc + '/' + 'Storage_%g.h5'%(time[k]), Storage)

                ########### flux

                write_hdf5(direc + '/' + 'Bflux_%g.h5'%(time[k]), Bflux)
                write_hdf5(direc + '/' + 'FPresflux_%g.h5'%(time[k]), FPressflux)
                write_hdf5(direc + '/' + 'PPresflux_%g.h5'%(time[k]), PPressflux)
                write_hdf5(direc + '/' + 'TauFflux_%g.h5'%(time[k]), TauFflux)
                write_hdf5(direc + '/' + 'TauPflux_%g.h5'%(time[k]), TauPflux)

                ############# bottom part

                write_hdf5(direc + '/' + 'B_l_%g.h5'%(time[k]), B_l)
                write_hdf5(direc + '/' + 'C_l_%g.h5'%(time[k]), C_l)
                write_hdf5(direc + '/' + 'TauFxy_l_%g.h5'%(time[k]), TauxyF_l)
                write_hdf5(direc + '/' + 'TauPxy_l_%g.h5'%(time[k]), TauxyP_l)
                write_hdf5(direc + '/' + 'dybdx_TauFxy_l_%g.h5'%(time[k]), dybdx_TauF_l)
                write_hdf5(direc + '/' + 'dybdx_TauPxy_l_%g.h5'%(time[k]), dybdx_TauP_l)
                write_hdf5(direc + '/' + 'fluidP_l_%g.h5'%(time[k]), FPres_l)
                write_hdf5(direc + '/' + 'particleP_l_%g.h5'%(time[k]), PPres_l)

                ############# top part

                write_hdf5(direc + '/' + 'B_h_%g.h5'%(time[k]), B_h)
                write_hdf5(direc + '/' + 'C_h_%g.h5'%(time[k]), C_h)
                write_hdf5(direc + '/' + 'TauFxy_h_%g.h5'%(time[k]), TauxyF_h)
                write_hdf5(direc + '/' + 'TauPxy_h_%g.h5'%(time[k]), TauxyP_h)
                write_hdf5(direc + '/' + 'dybdx_TauFxy_h_%g.h5'%(time[k]), dybdx_TauF_h)
                write_hdf5(direc + '/' + 'dybdx_TauPxy_h_%g.h5'%(time[k]), dybdx_TauP_h)
                write_hdf5(direc + '/' + 'fluidP_h_%g.h5'%(time[k]), FPres_l)
                write_hdf5(direc + '/' + 'particleP_h_%g.h5'%(time[k]), PPres_h)

        write_hdf5(direc + '/' + "xbed.h5", x)
        #write_hdf5("output_2D/y.h5", y)






