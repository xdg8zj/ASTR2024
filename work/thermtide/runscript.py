import math
import numpy as np
import os
import shutil
import sys

workdir = '/Users/aleynaloughran-pierce/Desktop/ASTR2024/work/thermtide'
athenadir = '/Users/aleynaloughran-pierce/Desktop/ASTR2024/athena-master'
pgendir = '/Users/aleynaloughran-pierce/Desktop/ASTR2024/athena-master/src/pgen'
pgenfile = 'thermtide.cpp'

def compile():
 shutil.copy(pgenfile,pgendir + '/' + pgenfile)
 os.chdir(athenadir)
 os.system('python3 configure.py --prob thermtide --coord spherical_polar --eos isothermal --ccmd /usr/bin/g++') 
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
 compile()
 run()
 #cleanpgenfile()

if __name__ == '__main__':
 main()
