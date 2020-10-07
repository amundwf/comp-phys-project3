import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os


### Optional: Run the C++ program to get an updated data file. 
# Compile and run the C++ files (this is exactly what is in the makefile):
os.system("echo compiling C++ codes...")
os.system("g++ -o main.out ../cpp_codes/main.cpp ../cpp_codes/utils.cpp -larmadillo")
os.system("echo executing...")
os.system("./main.out")
###

# Read the comma-separated data files (two columns, x and y):
directory = "../results/3a_earth_sun_system/"
filename = "earth_sun_euler.csv"
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


plot1, = plt.plot(xList, yList, label='forwardEuler')

plt.grid()
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.xlim(-1.5, 1.5)
plt.ylim(-1.5, 1.5)
plt.suptitle('Earth-Sun system, forward Euler')
plt.show()
