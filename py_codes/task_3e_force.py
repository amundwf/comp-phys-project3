import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os


runCppCode = True
if runCppCode == True:  
    # Compile and run the C++ files (this is exactly what is in the makefile):
    os.system("echo compiling C++ codes...")
    os.system("g++ -o main.out ../cpp_codes/main.cpp ../cpp_codes/utils.cpp ../cpp_codes/planet.cpp ../cpp_codes/solver.cpp -larmadillo")
    os.system("echo executing...")
    os.system("./main.out")


# Read the comma-separated data files (two columns, x and y):
directory = "../results/3e_force/"
betaFilePath = directory + "betaList.csv"

# Load the list of beta values:
betaValues = np.loadtxt(betaFilePath, skiprows=0, delimiter=",")
betaValues = pd.DataFrame(betaValues, columns=["beta"])
betaList = betaValues["beta"]
print(betaList)

betaListLength = len(betaList) # Number of beta values
# Retrieve all the results from files for the different beta values:
for i in range(betaListLength):
    # Match the file names in the results folder for 
    # task 3e:
    beta = betaList[i]
    betaStr = format(beta,'.2f')
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
    plot1, = plt.plot(xList, yList, label= "beta = " + betaStr, linewidth = 0.9)


# Plot the sun at the center of the solar system:
plt.plot(0, 0, 'r.', markersize=12, label = 'Sun')

plt.axis('equal')
plt.xlabel(r'$x$ (au)')
plt.ylabel(r'$y$ (au)')
lim = 2
plt.xlim(-lim, lim)
plt.ylim(-lim, lim)
plt.grid()
plt.legend()

plt.suptitle('The Sun-Earth system with varying force law')
plt.show()

 
