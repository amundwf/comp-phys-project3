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


# The results directory:
directory = "../results/3g_three_body/"


# Read the planet names from the planet_names.csv file:
'''
planetNamesData = np.loadtxt(directory + "planet_names.csv", skiprows=0, delimiter=",")
planetNamesData = pd.DataFrame(planetNamesData, columns=["planet_names"])
planetNames = planetNamesData["planet_names"]
'''
# Read in the file as a dataframe:
planetNamesDF = pd.read_csv(directory + "planet_names.csv", header=None, names=['planets'])
# Convert the dataframe to a list (of strings):
planetNames = planetNamesDF['planets'].to_list()

print(type(planetNames))
print("List of planet names:"); print(planetNames)
#planetNames = pd.read_csv(directory + "planet_names.csv")
#print("planetNames (python script):"); print(planetNames)

nPlanets = len(planetNames) # Number of planets, excluding the Sun

# Loop through all the planets:
for i in range(nPlanets):
    planetName = planetNames[i]
    # File name for the current planet:
    filename = planetName + ".csv"

    filePath = os.path.join(directory, filename) # The full file path.
    data = np.loadtxt(filePath, skiprows=1, delimiter=",")
    # Get the columns of data as lists:
    data = pd.DataFrame(data, columns=["t", "x", "y", "z", "vx", "vy", "vz"])
    tList = data["t"]
    xList = data["x"]
    yList = data["y"]
    zList = data["z"]
    vxList = data["vx"]
    vyList = data["vy"]
    vzList = data["vz"]

    # Plot the planet:
    print("Plotting orbit: " + planetName + " ...")
    plt.plot(xList, yList, label=planetName, linewidth = 0.8)
    

# Plot the sun at the center of the solar system:
plt.plot(0, 0, 'r.', markersize=12, label='Sun')

labelSize = 13
titleSize = 12

plt.xlabel(r'$x$ (au)', fontsize=labelSize)
plt.ylabel(r'$y$ (au)', fontsize=labelSize)
#plt.suptitle('Three-body problem, velocity Verlet', fontsize=titleSize)
#plt.xlim(-1.5, 1.5)
#plt.ylim(-1.5, 1.5)
plt.axis('equal')
plt.legend()
plt.grid()
plt.show()
 
