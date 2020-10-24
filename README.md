# comp-phys-project3 - Simulating the Solar System
Repository for project 3 in FYS4150 computational physics.

This is a project which simulated the motion of the planets in the Solar System. Both the Velocity Verlet and Euler's forward method area available in the code to be used. However, in the Velocity Verlet method is the only object oriented method. The object oriented code is composed of a Solver class to which planets can be added to and then run.  The Solver class does this by calling upon another class, Planet, to create planet object which it uses to update the position and velocity of the planets. Other functionaltiy exists for calculating the angular momentum, energy, percession of Mercury's perihelion and investigating the power law in Newton's equations. Each task is broken down in utils.cpp and each one can be run from main.cpp. The results are stored after the code runs in the results folder and can be plotted using the various python scripts. Some example plots are given in the results folder also. 

# In utils.cpp the various functions are as follows:

> Task 3a is for simulations of the Earth with Euler's method and Velocity Verlet without object orientation.

> Task 3b is for running simulations with object orientation of the Earth.

> Task 3d is for investigating different powers of r in the force equation calculation.

> Task 3f is for calculating escape velocities of planets.

> Task 3g is for simulating many bodies in the Solar system. 

> Task 3i is for calculating the perihelion shift of Mercury.

Although the code is broken into tasks in utils.cpp, the actual Solver and Planet class are very flexible for the problem you wish to experiment with. 

> Folders: cpp_codes, py_codes, results and report.
