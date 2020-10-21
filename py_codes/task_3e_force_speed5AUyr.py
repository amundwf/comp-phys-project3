import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os
import utils # utils.py

# This script performs the part of task 3e that asks for the elliptical
# orbit with initial speed 5 AU/yr and calculating the angular momentum.

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

# The mass of the Earth is required in order to calculate its angular momentum:
m_earth = 1

betaListLength = len(betaList) # Number of beta values.
# Save the lists of angular momentum (absolute values) in columns, one column
# for each beta value:
#angMomentumBetaArray = np.zeros((len(xList), betaListLength))

'''
# Retrieve all the results from files for the different beta values:
for i in range(betaListLength):
    # Match the file names in the results folder for 
    # task 3e:
    beta = betaList[i]
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

    # ###### FINISH THIS ANGULAR MOMENTUM STUFF ######
    print('Calculating angular momentum as a function of time...')
    angMomentumAbsList = utils.calculateAngularMomentumList(xList, yList, zList, vxList, vyList, vzList, m_e)
    print('Angular momentum calculation finished.')

    # Store the angular momentum for this beta value in column index i:
    angMomentumBetaArray[:,i] = angMomentumAbsList
    # plot(tList, angMomentumAbsList)

    # Plot each orbit, for all different values of beta:
    plot1, = plt.plot(xList, yList, label= "beta = " + betaStr)
'''

# Use the codeword to choose between plotting Earth's orbit (for all
# beta values) and plotting the angular momentum (for all beta values)
#codeword = 'orbit'
codeword = 'angularMomentum'

if codeword == 'orbit': # Plot the orbits:
    # Retrieve all the results from files for the different beta values:
    for i in range(betaListLength):
        # Match the file names in the results folder for 
        # task 3e:
        beta = betaList[i]
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
        plt.plot(xList, yList, label= "beta = " + betaStr)
        
    plt.axis('equal')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.xlim(-1.5, 1.5)
    plt.ylim(-1.5, 1.5)
    plt.suptitle('Sun-Earth system with varying force law, initial velocity 5 AU')
    # Plot the sun at the center of the solar system:
    plt.plot(0, 0, 'r.', markersize=20)

elif codeword == 'angularMomentum': # Plot angular momentum:
    for i in range(betaListLength): 
        beta = betaList[i]
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

        # Store the angular momentum for this beta value in column index i:
        angMomentumAbsList = utils.calculateAngularMomentumList(xList, yList, zList, vxList, vyList, vzList, m_earth)
        #angMomentumBetaArray[:,i] = angMomentumAbsList

        # Plot the angular momentum:
        plt.plot(tList, angMomentumAbsList)

    plt.xlabel(r'$t$')
    plt.ylabel(r'$L$')
    plt.suptitle('Angular momentum of the Sun-Earth system with varying force law, initial velocity 5 AU')

plt.grid()
plt.legend()
plt.show()


 
 
