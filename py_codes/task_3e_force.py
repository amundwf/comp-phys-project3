import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os


### Optional: Run the C++ program to get an updated data file. 
# Compile and run the C++ files (this is exactly what is in the makefile):
os.system("echo compiling C++ codes...")
os.system("g++ -o main.out ../cpp_codes/main.cpp ../cpp_codes/utils.cpp ../cpp_codes/planet.cpp ../cpp_codes/solver.cpp -larmadillo")
os.system("echo executing...")
os.system("./main.out")
###

# Read the comma-separated data files (two columns, x and y):
directory = "../results/3e_force/"

'''
betaListLength = 11 # Number of beta values
for i in range(1,betaListLength):

    # Fix the filename (using string formatting) to match the
    # filename in the results folder:   filename = "[something].csv"
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

    #plot1, = plt.plot(xList, yList)#, label='velocityVerlet')
'''

'''
# Plot the sun at the center of the solar system:
plt.plot(0, 0, 'r.', markersize=20)

plt.grid()
plt.axis('equal')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.xlim(-1.5, 1.5)
plt.ylim(-1.5, 1.5)
plt.suptitle('Sun-Earth system with varying force law (different values of beta)')
plt.show()
'''
 
