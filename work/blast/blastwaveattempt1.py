import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'/Users/aleynaloughran-pierce/Desktop/ASTR2024/athena-master/vis/python')
import athena_read
x,y,z,data = athena_read.vtk('/Users/aleynaloughran-pierce/Desktop/ASTR2024/athena-master/blastwave3/Blast.block0.out1.00007.vtk')
#print(data['rho'][:,1,1])
#radiusdensity = data['rho'][:,1,1]

plt.plot(x[0:len(x)-1],data['rho'][32,32,:]) #in order of phi, theta, r
plt.title("Density at all radaii")
plt.show()
