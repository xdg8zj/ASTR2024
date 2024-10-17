import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '../../athena-master/vis/python')
import athena_read
import plot_lines
import glob

#data = athena_read.hst('thermtide.hst')
#print('quadrapole moments {0}'.format(data['Qxx']))

def read_hst_file(filename):
  data = athena_read.hst(filename)
  return(data)
   #0 denotes the first argument in the 'data' section, while '1' denotes whatever follows

def plot_quadrupoles(filenames):
  for fname in filenames:
    data = read_hst_file(fname)
    time = data['time']
    dt = data['dt']
    mass = data['mass']
    mom_1 = data['1-mom']
    mom_2 = data['2-mom']
    mom_3 = data['3-mom']
    KE_1 = data['1-KE']
    KE_2 = data['2-KE']
    KE_3 = data['3-KE']
    tot_E = data['tot-E']
    Qxx = data['Qxx']
    Qxy = data['Qxy']
    Qxz = data['Qxz']
    Qyx = data['Qyx']
    Qyy = data['Qyy']
    Qyz = data['Qyz']
    Qzx = data['Qzx']
    Qzy = data['Qzy']
    Qzz = data['Qzz']
    plt.plot(time,Qxx, label = "Qxx", linestyle = "dotted")
#    plt.plot(time,Qxy, label = "Qxy", linestyle = "dotted")
#    plt.plot(time,Qxz, label = "Qxz", linestyle = "dotted")
#    plt.plot(time,Qyx, label = "Qyx", linestyle = "dotted")
    plt.plot(time,Qyy, label = "Qyy", linestyle = "dotted")
#    plt.plot(time,Qyz, label = "Qyz", linestyle = "dotted")
#    plt.plot(time,Qzx, label = "Qzx", linestyle = "dotted")
#    plt.plot(time,Qzy, label = "Qzy", linestyle = "dotted")
    plt.plot(time,Qzz, label = "Qzz", linestyle = "dotted")
  plt.yscale('log')
  plt.title('time of simulation vs log[Quadrupole moments]')
  plt.xlabel('time (s)')
  plt.ylabel('log [Quadrupole moments (kg*m^2)] ')
  plt.legend()
  plt.show()
  
def main():
  filenames = ['thermtide.hst']
  plot_quadrupoles(filenames)
  
if __name__ == '__main__':
  main()
   
