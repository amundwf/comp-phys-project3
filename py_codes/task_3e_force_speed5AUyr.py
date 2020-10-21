import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os
import utils # utils.py

# This script performs the part of task 3e that asks for the elliptical
# orbit with initial speed 5 AU/yr and calculating the angular momentum.
# IMPORTANT: Make sure that the initial speed in task_3e_force() in utils.cpp is
# set to 5 AU/yr.
# Choose which beta values to run for in task_3e_force() in utils.cpp.

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
betaValues = np.loadtxt(betaFilePath, skiprows=0, delimiter=",")
betaValues = pd.DataFrame(betaValues, columns=["beta"])
betaList = betaValues["beta"]

# The masses of the Earth and the Sun are required in order to calculate the 
# total energy and angular momentum of the Earth:
m_earth = 3.002463426e-6 # In solar masses.
m_sun = 1

betaListLength = len(betaList) # Number of beta values.
# Save the lists of angular momentum (absolute values) in columns, one column
# for each beta value:
#angMomentumBetaArray = np.zeros((len(xList), betaListLength))

# Use the codeword to choose between plotting Earth's orbit (for all
# beta values) and plotting the angular momentum (for all beta values)
#codeword = 'orbits'
#codeword = 'totalEnergy'
codeword = 'angularMomentum'

if codeword == 'orbits': # Plot the orbits:
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
        plt.plot(xList, yList, label= "beta = " + betaStr, linewidth = 0.9)   
    plt.axis('equal')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.xlim(-1.5, 1.5)
    plt.ylim(-1.5, 1.5)
    plt.suptitle('Sun-Earth system with varying force law, initial velocity 5 AU')
    # Plot the Sun at the center of the solar system:
    plt.plot(0, 0, 'r.', markersize=15, label = 'Sun')

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

        # Calculate the total energy for all points in time:
        angMomentumAbsList = utils.calculateAngularMomentumList(xList, yList, zList, \
            vxList, vyList, vzList, m_earth)
        # Plot the angular momentum:
        plt.plot(tList, angMomentumAbsList, label= "beta = " + betaStr, linewidth = 0.8)
    # y axis limits: 
    minimum = min(angMomentumAbsList)
    maximum = max(angMomentumAbsList)
    diff = maximum - minimum
    avg = (maximum + minimum)/2 # 'average'
    plt.ylim(avg - diff*1.4, avg + diff*1.4)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$L$')
    plt.suptitle('Angular momentum of the Sun-Earth system with \nvarying force law, initial velocity 5 AU/yr')

elif codeword == 'totalEnergy':
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

        # Calculate the total energy for all points in time:
        totalEnergyList = utils.calculateTotalEnergyList(xList, yList, zList, \
            vxList, vyList, vzList, m_earth, m_sun)
        # Plot the total energy:
        plt.plot(tList, totalEnergyList, label= "beta = " + betaStr)
    # y axis limits:
    minimum = min(totalEnergyList)
    maximum = max(totalEnergyList)
    diff = maximum - minimum
    avg = (maximum + minimum)/2 # 'average'
    print("avg value:"); print(avg)
    plt.ylim(avg - diff, avg + diff)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')
    plt.suptitle('Total energy of the Sun-Earth system with varying\n force law, initial velocity 5 AU/yr')

plt.grid()
plt.legend()
plt.show()


 
 
