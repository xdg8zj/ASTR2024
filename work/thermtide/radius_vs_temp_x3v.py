#read vtk file into python:

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '../../athena-master/vis/python')
import athena_read
import plot_spherical
import glob


def read_vtk_file(filename):
    x,y,z,data = athena_read.vtk(filename)
 #print("headers=",data.keys())
    return x,y,z,data


def plot_rad_temp(filenames):
    for fname in filenames:
        x,y,z,data = read_vtk_file(fname)
        press = data['press']
        rho = data['rho']
        mu = 2.3
        amu = 1.67262192595e-24
        kboltz = 1.3807e-16
        ndof = 5.0
        
        nk = len(z) -1
        nj = len(y)-1
        ni = len (x)-1
        xc=0.5*(x[0:ni]+np.roll(x,-1)[0:ni])
        
        k = 0
        j = 0
        

        
        
        for k in range(nk):
            for j in range(nj):
                temp = mu*amu*press[k,j,:]/rho[k,j,:]/kboltz
            plt.scatter(xc,temp)
            
        plt.xlim(7.0e9,7.1e9)
        plt.legend()
        plt.show()


#        for k in z:
#            for j in y:
#                for i in x:
#                    mu = 2.3
#                    amu = 1.67262192595e-24
#                    kboltz = 1.3807e-16
#                    ndof = 5.0
#                    
#                    temp[k,j,i] =(press[k,j,i])/rho[k,j,i]

##                    
#                    
#        print(temp)
#        temp_rad = temp[0,0,:]
#        plt.plot(rad,temp_rad, label = "Radius temp", linestyle = "dotted")
#        plt.yscale('log')
#        plt.title('radius (m) vs temp (K)')
#        plt.xlabel('radius (m)')
#        plt.ylabel('log[Temp (K)]')
#        plt.show()
#        
def main():
    filenames = ['thermtide.block0.out1.02000.vtk']
    plot_rad_temp(filenames)
    
    
if __name__ == '__main__':
    main()
        
   
   

