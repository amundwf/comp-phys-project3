import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os
import time
import utils # utils.py

# This script performs the part of task 3e that asks for the elliptical
# orbit with initial speed 5 AU/yr and calculating the angular momentum.
# IMPORTANT: Make sure that the initial speed in task_3e_force() in utils.cpp is
# set to 5 AU/yr.
# Choose which beta values to run for in task_3e_force() in utils.cpp.

# Choose whether to run the C++ code. If the results files already exist and
# only plot tweaking was needed, you can set runCppCode to False to skip
# running the C++ code. Set runCppCode to True if you want to run the C++
# code.
runCppCode = True
if runCppCode == True:
    # Compile and run the C++ files (this is exactly what is in the makefile):
    os.system("echo compiling C++ codes...")
    tic = time.perf_counter()
    os.system("g++ -o main.out ../cpp_codes/main.cpp ../cpp_codes/utils.cpp ../cpp_codes/planet.cpp ../cpp_codes/solver.cpp -larmadillo")
    toc = time.perf_counter()
    print(f"C++ compilation took {toc - tic:0.3f} seconds")

    os.system("echo executing...")
    # Time the C++ execution time:
    tic = time.perf_counter()
    os.system("./main.out") # Execute the C++ code
    toc = time.perf_counter()
    print(f"C++ execution took {toc - tic:0.3f} seconds")
    


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
#codeword = 'totalEnergyAngMomentum'
codeword = 'totalEnergy'
#codeword = 'angularMomentum'
#codeword = 'distance' # Plot the distance from the Sun as a function of time.

print('Plotting...')
tic = time.perf_counter() # Time the plotting
labelSize = 13
titleSize = 12
if codeword == 'orbits': # Plot the orbits:
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
        tList = data["t"]
        xList = data["x"]
        yList = data["y"]
        zList = data["z"]
        vxList = data["vx"]
        vyList = data["vy"]
        vzList = data["vz"]

        # Plot each orbit, for all different values of beta:
        plt.plot(xList, yList, label= "beta = " + betaStr, linewidth = 1.1)   
    plt.axis('equal')
    plt.xlabel(r'$x$', fontsize=labelSize)
    plt.ylabel(r'$y$', fontsize=labelSize)
    plt.xlim(-1.5, 1.5)
    plt.ylim(-1.5, 1.5)
    plt.suptitle('Sun-Earth system with varying force law, initial velocity 5 au/yr')
    # Plot the Sun at the center of the solar system:
    plt.plot(0, 0, 'r.', markersize=15, label = 'Sun')
    plt.legend()

elif codeword == 'totalEnergyAngMomentum': # Plot angular momentum:
    for i in range(betaListLength): 
        beta = betaList[i]
        betaStr = format(beta,'.2f')
        filename = '3e_force_beta' + betaStr + '.csv'
        filePath = os.path.join(directory, filename) # The full file path.

        # Load the data files:
        data = np.loadtxt(filePath, skiprows=1, delimiter=",")
        # Get the columns as lists:
        data = pd.DataFrame(data, columns=["t", "x", "y", "z", "vx", "vy", "vz"])
        tList = data["t"]
        xList = data["x"]
        yList = data["y"]
        zList = data["z"]
        vxList = data["vx"]
        vyList = data["vy"]
        vzList = data["vz"]

        # Calculate the total energy for all points in time:
        angMomentumAbsList = utils.calculateAngularMomentumList(xList, yList, zList, \
            vxList, vyList, vzList, m_earth)
        L0 = angMomentumAbsList[0] # Initial angular momentum
        angMomentumNormalized = angMomentumAbsList/L0
        # Calculate the total energy for all points in time:
        totalEnergyList = utils.calculateTotalEnergyList(xList, yList, zList, \
            vxList, vyList, vzList, m_earth, m_sun)
        E0 = totalEnergyList[0] # Initial total energy
        totalEnergyNormalized = totalEnergyList/E0

        # Plot:
        plt.plot(tList, totalEnergyNormalized, label= "E, beta = " + betaStr, linewidth = 0.8)
        plt.plot(tList, angMomentumNormalized, label= "L, beta = " + betaStr, linewidth = 0.8)
    # y axis limits: 
    #minimum = min(totalEnergyNormalized)
    #maximum = max(totalEnergyNormalized)
    #diff = maximum - minimum
    #avg = (maximum + minimum)/2 # 'average'
    #plt.ylim(avg - diff*1.4, avg + diff*1.4)
    plt.xlabel(r'$t$ (yr)', fontsize=labelSize)
    plt.ylabel(r'$E/E0, L/L0$', fontsize=labelSize)
    plt.suptitle('Total energy and angular momentum of the Sun-Earth system with \nvarying force law, initial velocity 5 au/yr')
    plt.legend()

elif codeword == 'totalEnergy':
    # Plot beta=2 last to get it in the forefront:
    #for i in range(betaListLength-1,-1,-1):
    for i in range(betaListLength):
        beta = betaList[i]
        betaStr = format(beta,'.2f')
        filename = '3e_force_beta' + betaStr + '.csv'
        filePath = os.path.join(directory, filename) # The full file path.

        # Load the data files:
        data = np.loadtxt(filePath, skiprows=1, delimiter=",")
        # Get the columns as lists:
        data = pd.DataFrame(data, columns=["t", "x", "y", "z", "vx", "vy", "vz"])
        tList = data["t"]
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
        lw = 1.3 # line width
        if beta==2.0:
            lw = 2.5
        plt.plot(tList, totalEnergyList, label= "beta = " + betaStr, linewidth = lw)
        #plt.semilogy(tList, totalEnergyList, label= "beta = " + betaStr)
    # y axis limits:
    minimum = min(totalEnergyList)
    maximum = max(totalEnergyList)
    diff = maximum - minimum
    avg = (maximum + minimum)/2 # 'average'
    #print("avg value:"); print(avg)
    #plt.ylim(avg - diff, avg + diff)
    plt.xlabel(r'$t$ (yr)', fontsize=labelSize)
    plt.ylabel(r'$E$', fontsize=labelSize)
    plt.suptitle('Total energy of the Sun-Earth system with varying\n force law, initial velocity 5 au/yr')
    plt.legend()

elif codeword == 'angularMomentum': # Plot angular momentum:
    for i in range(betaListLength): 
        beta = betaList[i]
        betaStr = format(beta,'.2f')
        filename = '3e_force_beta' + betaStr + '.csv'
        filePath = os.path.join(directory, filename) # The full file path.

        # Load the data files:
        data = np.loadtxt(filePath, skiprows=1, delimiter=",")
        # Get the columns as lists:
        data = pd.DataFrame(data, columns=["t", "x", "y", "z", "vx", "vy", "vz"])
        tList = data["t"]
        xList = data["x"]
        yList = data["y"]
        zList = data["z"]
        vxList = data["vx"]
        vyList = data["vy"]
        vzList = data["vz"]

        # Calculate the total energy for all points in time:
        tic = time.perf_counter() # Time the calculations
        angMomentumAbsList = utils.calculateAngularMomentumList(xList, yList, zList, \
            vxList, vyList, vzList, m_earth)
        L0 = angMomentumAbsList[0] # Initial angular momentum
        angMomentumNormalized = angMomentumAbsList/L0
        toc = time.perf_counter()
        print(f"Python calculation took {toc - tic:0.3f} seconds")

        # Plot:
        plt.plot(tList, angMomentumAbsList, label= "beta = " + betaStr, linewidth = 1.2)
        #plt.semilogy(tList, angMomentumAbsList, label= "beta = " + betaStr, linewidth = 1.0)
        toc = time.perf_counter()
        print(f"Plotting took {toc - tic:0.3f} seconds")
    # y axis limits:
    minimum = min(angMomentumAbsList)
    maximum = max(angMomentumAbsList)
    diff = maximum - minimum
    avg = (maximum + minimum)/2 # 'average'
    limitMultiplier = 2
    plt.ylim(avg - diff*limitMultiplier, avg + diff*limitMultiplier)
    plt.xlabel(r'$t$ (yr)', fontsize=labelSize)
    plt.ylabel(r'$L$', fontsize=labelSize)
    plt.suptitle('Angular momentum of the Sun-Earth system with \nvarying force law, initial velocity 5 au/yr')
    plt.legend()

elif codeword == 'distance':
    for i in range(betaListLength): 
        beta = betaList[i]
        betaStr = format(beta,'.2f')
        filename = '3e_force_beta' + betaStr + '.csv'
        filePath = os.path.join(directory, filename) # The full file path.

        # Load the data files:
        data = np.loadtxt(filePath, skiprows=1, delimiter=",")
        # Get the columns as lists:
        data = pd.DataFrame(data, columns=["t", "x", "y", "z", "vx", "vy", "vz"])
        tList = data["t"]
        xList = data["x"]
        yList = data["y"]
        zList = data["z"]
        vxList = data["vx"]
        vyList = data["vy"]
        vzList = data["vz"]

        # List of the absolute values of the 3D positions contained in the position
        # vectors from xList, yList and zList.
        distanceList = utils.calculateDistanceList(xList, yList, zList)

        # Plot the distance in time:
        plt.semilogy(tList, distanceList, label = "beta = " + betaStr, linewidth = 0.9)
    # y axis limits:
    minimum = min(distanceList)
    maximum = max(distanceList)
    diff = maximum - minimum
    avg = (maximum + minimum)/2 # 'average'
    print("avg value:"); print(avg)
    plt.ylim(avg - diff, avg + diff)
    
    plt.xlabel(r'$t$ (yr)', fontsize=labelSize)
    plt.ylabel(r'$r$ (au)', fontsize=labelSize)
    plt.suptitle('The distance between the Earth and the Sun with varying\n\
         force law, initial velocity 5 au/yr', fontsize=titleSize)
    plt.legend(loc = 'upper right')

plt.grid()
toc = time.perf_counter()
print(f"Plotting took {toc - tic:0.3f} seconds")

plt.show()


 
 
