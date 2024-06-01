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




#class Budget_Integrand:

#    def __init__(self, xbed, Yb):

#        self.xbed = xbed
#        self.Yb   = Yb
        #self.alpha = alpha
        #self.ybedmethod = ybedmethod


#    def bed_height(self):

#        bed_l   = np.zeros(len(self.xbed), dtype = int)
#        bed_h   = np.zeros(len(self.xbed), dtype = int)

#        ybed_l   = np.zeros(len(self.xbed))
#        ybed_h   = np.zeros(len(self.xbed))
        
#        _, ny, _ = self.Yb.shape
#        alphamin = 1e-3
#        per = 0.9

       # xmax = np.max(self.xbed)
       # ymax = np.max(self.Yb[0, :, 0])

       # (n1d, ny, n2d) = np.shape(self.Yb)

       # for i in range(n1d):

                        
        #    bed_l[i] = 0
        #    bed_h[i] = ny-1


         #   ybed_l[i] = self.Yb[i, bed_l[i], 0]
          #  ybed_h[i] = self.Yb[i, bed_h[i], 0]

       # return bed_l, bed_h, ybed_l, ybed_h

    
   # def x_deri_1D(self, A, x):
        
    #    deri_x = np.zeros(A.shape)

     #   for i in range(len(x)):
      #      if i==0:
       #         deri_x[i] = (A[i+1] - A[i]) / (x[i+1] - x[i])

        #    elif 0<i<len(x)-1:
         #       deri_x[i] = (A[i+1] - A[i-1]) / (x[i+1] - x[i-1])

          #  else:
           #     deri_x[i] = (A[i] - A[i-1]) / (x[i] - x[i-1])
       # return deri_x


   # def convective_flux(self, A, dybdx, bed):

    #    nx, ny, nz = A.shape

     #   B = np.zeros(nx)

      #  for i in range(len(bed)):
       #     B[i] = A[i, bed[i],0]* dybdx[i]
        #return B

   # def Integrate_flux(self, A, Yb, bed_l, bed_h):

    #    nx, ny, nz = A.shape

     #   InteValue = np.zeros(nx)
      #  for i in range(nx):
       #     integrand = A[i, bed_l[i]:bed_h[i], 0]
        #    Ybinte    = Yb[i, bed_l[i]: bed_h[i], 0]
         #   InteValue[i] = np.trapz(integrand, Ybinte)
       # return InteValue

   # def Surfacevalue(self, A, bed):

    #    nx, ny, nz = A.shape
     #   B = np.zeros(nx)
        
      #  for i in range(len(bed)):
       #     B[i] = A[i, bed[i], 0]
       # return B

   # def SurfacevalueMul(self, A,B, bed):
 #       
  #      nx, ny, nz = A.shape
 #       C = np.zeros(nx)
#
  #      for i in range(len(bed)):
 #           C[i] = A[i, bed[i], 0]* B[i]
#        return C



class Budget_Integrand:

    def __init__(self, xbed, Yb, alpha=None, per=None, lower=None, upper = None):
        self.xbed = xbed
        self.Yb   = Yb
        self.alpha = alpha
        self.lower = lower
        self.upper = upper
        self.per   = per

    def bed_height(self):
        
        xmax = np.max(self.xbed)
        ymax = np.max(self.Yb[0, :, 0])
        (nx, ny, nz) = np.shape(self.Yb)
        
        nx, nz = np.shape(self.xbed)
        
        bed_l   = np.zeros((nx, nz), dtype = int)
        bed_h   = np.zeros((nx, nz), dtype = int)
        ybed_l   = np.zeros((nx, nz))
        ybed_h   = np.zeros((nx, nz))
        alphamin = 1e-3

        alpha_shape = np.amax(self.alpha)

        for i in range(nx):

            for j in range(nz):


                if alpha_shape!=0:
    
                    if self.lower!=None and self.upper!=None:
                        alphamax = np.round(self.per*np.amax(self.alpha[i, :, j]),3)

                        if self.xbed[i]>0 and self.xbed[i]<=xmax - 0.5*ymax:
                            if np.size(np.where(np.round(self.alpha[i, :, j],3)>=alphamax)[0])==0:
                                bed_l[i, j] = 0
                                bed_h[i, j] = 0
                            else:
                                index = np.where(np.round(self.alpha[i, :, j], 3) > alphamax)
                                if len(index[0])>0:
                                    bed_l[i, j] = index[0][-1]

                                index_h = np.where(np.round(self.alpha[i, :, j],3)> alphamin)

                                if len(index_h[0])>0:
                                    bed_h[i, j] = index_h[0][-1]
                    else:
                        bed_l[i, j] = 0
                        bed_h[i, j] = ny - 1
                else:
                    print ("Manohar is in without alpha loop")
                    bed_l[i, j] = 0
                    bed_h[i, j] = ny - 1

                ybed_l[i, j] = self.Yb[i, bed_l[i, j], j]
                ybed_h[i, j] = self.Yb[i, bed_h[i, j], j]

        return bed_l, bed_h, ybed_l, ybed_h

    def x_deri(self, A, x):

        #print ("x_deri = \t", x.shape)

        deri_x = np.zeros(A.shape)

        x = x[:, 0]
        #print ("x_deri = \t", x.shape)

        for i in range(len(x)):
                
            if i==0:
                deri_x[i, :] = (A[i+1, :] - A[i, :]) / (x[i+1] - x[i])

            elif 0<i<len(x)-1:
                deri_x[i, :] = (A[i+1, :] - A[i-1, :]) / (x[i+1] - x[i-1])

            else:
                deri_x[i, :] = (A[i, :] - A[i-1, :]) / (x[i] - x[i-1])

        return deri_x

    def z_deri_2D(self, A, z):

        deri_z = np.zeros(A.shape)
        
        z = z[0, :]

        for i in range(len(z)):

            if i==0:
                deri_z[:, i] = (A[:, i+1] - A[:, i]) / (z[i+1] - z[i])

            elif 0<i<len(x)-1:
                deri_z[:, i] = (A[:, i+1] - A[:, i-1]) / (z[i+1] - z[i-1])

            else:
                deri_z[:, i] = (A[:, i] - A[:, i-1]) / (z[i] - z[i-1])
        return deri_z


    def Integrate_flux(self, A, Yb, bed_l, bed_h):
        
        nx, ny, nz = A.shape
        InteValue = np.zeros((nx, nz))

        for i in range(nx):
            for j in range(nz):
                integrand = A[i, bed_l[i, j]:bed_h[i, j], j]
                Ybinte    = Yb[i, bed_l[i, j]: bed_h[i, j], j]
                InteValue[i, j] = np.trapz(integrand, Ybinte)
        return InteValue


    def Surfacevalue(self, A, bed):

        nx, ny, nz = A.shape
        B = np.zeros((nx, nz))

        for i in range(nx):
            for j in range(nz):
                B[i] = A[i, bed[i, j], j]
        return B

    def SurfacevalueMul(self, A,B, bed):

        nx, ny, nz = A.shape
        C = np.zeros((nx, nz))

        for i in range(nx):
            for j in range(nz):
                C[i, j] = A[i, bed[i, j], j]* B[i, j]
        return C





























