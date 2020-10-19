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
betaFilePath = directory + "betaList.csv"

# Load the list of beta values:
betaValues = np.loadtxt(betaFilePath, skiprows=1, delimiter=",")
betaValues = pd.DataFrame(betaValues, columns=["beta"])
betaList = betaValues["beta"]
print(betaList)

betaListLength = len(betaList) # Number of beta values
# Retrieve all the results from files for the different beta values:
for beta in betaList:
    # Match the file names in the results folder for 
    # task 3e:
    betaStr = format(beta,'.1f')
    filename = '3e_force_beta' + betaStr + '.csv'
    filePath = os.path.join(directory, filename) # The full file path.

    # Load the data files:
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

    # Plot each orbit, for all different values of beta:
    plot1, = plt.plot(xList, yList, label= "beta = " + betaStr)


# Plot the sun at the center of the solar system:
plt.plot(0, 0, 'r.', markersize=20)

plt.grid()
plt.legend()
plt.axis('equal')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.xlim(-1.5, 1.5)
plt.ylim(-1.5, 1.5)

plt.suptitle('Sun-Earth system with varying force law (different values of beta)')
plt.show()

 
