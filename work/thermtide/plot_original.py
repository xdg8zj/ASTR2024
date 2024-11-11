# read vtk file into python: https://github.com/PrincetonUniversity/athena/wiki/Reading-Data-into-Python

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

def plot_den_x(filenames):
 for fname in filenames:
   x,y,z,data = read_vtk_file(fname)
   rho=data['rho']
   x1rho = data['x1rho']
    
   n=len(x)-1
   xc=0.5*(x[0:n]+np.roll(x,-1)[0:n])
   k=0
   j=0
   plt.plot(xc,rho[k,j,:], label = "simulation", color = "blue", linestyle = "dotted")
 plt.yscale('log')
# plt.title('radius and density (actual data)')
 plt.xlabel('radius (cm)')
 plt.ylabel('density (g/cm^3)')
 plt.legend()
 plt.show()
 
 
    
    

def main():
 #filenames = glob.glob('*.vtk')
 filenames = ['thermtide.block0.out1.01000.vtk']
 #plot_den_x(filenames)
 #plt.close()
 for fname in filenames:
    x,y,z,data = read_vtk_file(fname)
    rho = data['rho']
    x1rho = data['x1rho']
    
 rho_base = rho[0,0,0]
 r_base = x[0]
 print(x1rho)
# print(rho_base, r_base)
#
# 
# #plotting another plot for comparison
# 
# radius_compare_lists = x
# print(radius_compare_lists)
#
# GM = (6.67e-8)*(1.89813e30)
# gamma_compare = 1.4
# P_const = 10e6
# rho_const =(10e-4)
# 
# #for K critical:
# K_crit = GM*(gamma_compare-1)/(gamma_compare*r_base*(rho_base**(gamma_compare-1)))
## C = 0
# print("K_crit:", K_crit)
# print("gamma:", gamma_compare)
## quit()
## K_bottom = rho_const**(gamma_compare)
## K_compare = P_const/K_bottom
## K_larger_than_crit = 2*K_crit
## K_equals_crit = K_crit
## K_smaller_than_crit = K_crit/5
#
# K_ratio = 0.0147
#
# radius_ratio = r_base/radius_compare_lists
# print("radius ratio:", radius_ratio)
## C = (1/r_base)*((K_smaller_than_crit/K_crit)-1)
# rho_compare = rho_base*((1/K_ratio)*(radius_ratio+K_ratio-1))**(1/(gamma_compare-1))
## rho_compare = (((1/radius_compare_lists)+C)*GM*(gamma_compare-1)/(K_smaller_than_crit*gamma_compare))**(1/(gamma_compare-1))
# print("New rho:", rho_compare)
# plot_compare = plt.plot(radius_compare_lists, rho_compare, label = f"formula plot \nK crit = {K_crit} \nK/Kcrit ratio = {K_ratio}", color = "red")
# plt.plot(plot_compare,plot_den_x(filenames))
# plt.yscale('log')
# plt.xlabel('radius (cm)')
# plt.ylabel('density (g/cm^3)')
# plt.title('desnity function graph')
## plt.text(7.1e9,10e-4, f"P = {P_const} \n rho = {rho_bottom} \n K = {K_compare}")
# plt.legend()
# plt.show()
## plt.savefig(f"k_kcrit_{K_ratio}.png")
# 
# 
# 
# 
#
#if __name__ == '__main__':
# main()
