import math
import numpy as np
import os
import shutil
import sys

workdir = '/Users/aleynaloughran-pierce/Desktop/ASTR2024/work/thermtide'
athenadir = '/Users/aleynaloughran-pierce/Desktop/ASTR2024/athena-master'
pgendir = '/Users/aleynaloughran-pierce/Desktop/ASTR2024/athena-master/src/pgen'
pgenfile = 'thermtide.cpp'

def configure():
 #shutil.copy(pgenfile,pgendir + '/' + pgenfile)
 os.chdir(athenadir)
 # serial
 os.system('python3 configure.py --prob=thermtide --coord=spherical_polar --ccmd=/usr/bin/g++') 
 # parallel with open-mpi and also  hdf5. this compiles! still need to test.
 #os.system('python3 configure.py -mpi -hdf5 --prob=thermtide --coord=spherical_polar --mpiccmd=/opt/homebrew/Cellar/open-mpi/5.0.5/bin/mpic++ --hdf5_path=/opt/homebrew/Cellar/hdf5-mpi/1.14.4.3') 
 # parallel with open-mpi, no hdf5. this compiles! still need to test.
 #os.system('python3 configure.py -mpi --prob=thermtide --coord=spherical_polar --mpiccmd=/opt/homebrew/Cellar/open-mpi/5.0.5/bin/mpic++')
 os.chdir(workdir)

def compile():
 shutil.copy(pgenfile,pgendir + '/' + pgenfile)
 os.chdir(athenadir)
 #os.system('python3 configure.py --prob=thermtide --coord=spherical_polar --ccmd=/usr/bin/g++') 
 #os.system('python3 configure.py -mpi -hdf5 --prob=thermtide --coord=spherical_polar --mpiccmd=/opt/homebrew/bin/mpicc --hdf5_path=/opt/homebrew/Cellar/hdf5-mpi/1.14.3') 
 #os.system('python3 configure.py -mpi --prob=thermtide --coord=spherical_polar --mpiccmd=/opt/homebrew/bin/mpicc') 
 os.system('make clean')
 os.system('make')
 os.chdir(workdir)

def run():
 os.system(athenadir + '/bin/athena' + ' -i athinput.thermtide')

def cleandata():
 os.system('rm -rf *.vtk')

def cleanpgenfile():
 os.system('rm ' + pgendir + '/' + pgenfile)

def main():

 configure()
 compile()
 run()
 #cleanpgenfile()

if __name__ == '__main__':
 main()
