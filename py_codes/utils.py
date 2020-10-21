# This file contains useful functions to be used in other .py files.

import numpy as np

def calculateAngularMomentumList(xList, yList, zList, vxList, vyList, vzList, mass):
    # Given data lists containing positions and velocities (x,y,z)
    # as functions of time, this function returns a vector containing the
    # corresponding angular momentum absolute values at each timestep.
    # Also, the mass of the object (planet) is a required input.
    length = len(xList)
    angMomentumAbsList = np.zeros(length)

    # Calculate the angular momentums. In order to avoid using too much memory
    # by creating large matrices, only one-dimensional lists are used.
    for i in range(length):
        position = np.array([xList[i], yList[i], zList[i]])
        velocity = np.array([vxList[i], vyList[i], vzList[i]])
        angMomentum = mass * np.cross(position,velocity)
        # Get the absolute value and store it in the vector:
        angMomentumAbsList[i] = np.linalg.norm(angMomentum)

    return angMomentumAbsList

def calculateTotalEnergyList(xList, yList, zList, vxList, vyList, vzList, massThis, massOther):
    # This function calculates the total energy (potential + kinetic) of an object
    # in orbit around another object. This object has mass massThis (the Earth), and the other
    # object has mass massOther (the Sun).
    length = len(xList)
    totalEnergyList = np.zeros(length)

    # Calculate the total energies. In order to avoid using too much memory
    # by creating large matrices, only one-dimensional lists are used.
    G = 4*(np.pi**2) # The gravitational constant in units AU, years and solar masses.
    for i in range(length):
        # Distance vector between the two objects:
        distanceVec = np.array([xList[i], yList[i], zList[i]])
        # Velocity of this object:
        velocity = np.array([vxList[i], vyList[i], vzList[i]])
        # Get the absolute values (norms) of the vectors:
        r = np.linalg.norm(distanceVec, ord=2)
        v = np.linalg.norm(velocity, ord=2)

        # Calculate the potential (U), kinetic (K) and total (E) energies:
        U = -(G*massThis*massOther)/(r**2)
        K = 0.5*massThis*(v**2)
        E = U + K
        totalEnergyList[i] = E

    return totalEnergyList