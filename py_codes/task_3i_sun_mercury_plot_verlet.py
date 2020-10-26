import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os


### Optional: Run the C++ program to get an updated data file. 
# Compile and run the C++ files (this is exactly what is in the makefile):
'''
os.system("echo compiling C++ codes...")
os.system("g++ -o main.out ../cpp_codes/main.cpp ../cpp_codes/utils.cpp -larmadillo")
os.system("echo executing...")
os.system("./main.out")
###
'''
# Read the comma-separated data files (two columns, x and y):
directory = "../results/3i_mercury_precession/"

############################################################################
filename = "mercury_sun_verlet.csv"
#filePath = fullfile(directory, filename)
filePath = os.path.join(directory, filename) # The full file path.
data = np.loadtxt(filePath, skiprows=1, delimiter=",")

# Get the columns as lists:
data = pd.DataFrame(data, columns=["t", "x", "y", "z", "vx", "vy", "vz"])

tList = data["t"] # Same for all values of omega_r
xList = data["x"]
yList = data["y"]
zList = data["z"]
vxList = data["vx"]
vyList = data["vy"]
vzList = data["vz"]

plot1, = plt.plot(xList, yList, label='velocityVerlet')

plt.grid()
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.xlim(-0.5, 0.5)
plt.ylim(-0.5, 0.5)
plt.suptitle('Mercury-Sun system with relativistic correction')
plt.savefig(directory + "/mercury-sun.pdf", dpi = "figure")
plt.show()

#############################################################################

filename = "perihelion.csv"
filePath = os.path.join(directory, filename) 
data = np.loadtxt(filePath, skiprows=1, delimiter=",")

# Get the columns as lists:
data = pd.DataFrame(data, columns=["x", "y", "z"])
data = data.loc[~(data==0).all(axis=1)]

############################################################################
fig = plt.figure()
plt.scatter(data["x"], data["y"], marker = ".", linewidth = 0.1)
plt.grid()
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.xlim(data["x"].min(), data["x"].max() )
plt.ylim(data["y"].min(), data["y"].max() )
plt.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
plt.suptitle('Position of the perihelion over a period of a century')
plt.savefig(directory + "/perihelion.pdf", dpi = "figure")
plt.show()

#############################################################################
fig = plt.figure()

theta = np.arctan(data["y"]/data["x"])*360*3600/(2*np.pi)
plt.grid()
plt.ylabel("perihelion angle /arc seconds")
plt.xlabel("revolutions")

x = np.arange(len(theta))
plt.scatter(x, theta, marker = ".", linewidth = 0.1)
plt.savefig(directory + "/angle_theta.pdf", dpi = "figure")
print(data.shape)
print(theta)
