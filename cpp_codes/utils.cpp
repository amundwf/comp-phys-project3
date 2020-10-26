#include <iostream>
#include <cmath>
#include <armadillo>
#include <sstream>  // std::stringstream objects
#include <iomanip>  // std::setprecision?
#include "utils.hpp"
#include "solver.hpp" 
#include "planet.hpp"

using namespace std;
using namespace arma;

// This file contains various functions, for running different parts of the project.

vec gForceVector(double G, double mass1, double mass2, vec pos1, vec pos2){
    // This returns the gravitational force *on* object 1 *from* object 2. Object 1
    // is pulled towards object 2.
    // Example: In sun-earth problem for the Earth orbit, the Earth is object 1 and
    // the Sun is object 2. 
    // G: gravitational constant.
    double r = norm(pos2 - pos1); // Relative distance between the objects.

    double forceStrength = (G*mass1*mass2)/(r*r);   // Newton's gravitational law.

    vec forceDirection = (pos2-pos1)/norm(pos2-pos1);   // This vector points *from*
    // object 1 and *towards* object 2, meaning that object 1 is influenced by object 2.
    vec forceVector = forceStrength * forceDirection;
    //cout << forceVector.t() << endl;
    return forceVector;
}

vec gForceVectorPlanet(Planet planet1, Planet planet2, double G){
    // This returns the gravitational force *on* planet 1 *from* planet 2. Object 1
    // is pulled towards object 2.
    // Example: In sun-earth problem for the Earth orbit, the Earth is object 1 and
    // the Sun is object 2.
    // G: gravitational constant.
    // Units: length: au, mass: solar masses, time: years.
    vec pos1 = planet1.position;
    double mass1 = planet1.mass;
    vec pos2 = planet2.position;
    double mass2 = planet2.mass;

    double r = norm(pos2 - pos1); // Relative distance between the objects.

    double forceStrength = (G*mass1*mass2)/(r*r);   // Newton's gravitational law.

    vec forceDirection = (pos2-pos1)/norm(pos2-pos1);   // This vector points *from*
    // object 1 and *towards* object 2, meaning that object 1 is influenced by object 2.
    vec forceVector = forceStrength * forceDirection;
    //cout << forceVector.t() << endl;
    return forceVector;
}

vec gForcePlanetBeta(Planet planet1, Planet planet2, double beta, double G){
    // Testing other forms of the gravitational force from the inverse-square
    // law.
    // This returns the gravitational force *on* planet 1 *from* planet 2. Object 1
    // is pulled towards object 2.
    // G: gravitational constant.
    // Units: length: au, mass: solar masses, time: years.
    vec pos1 = planet1.position;
    double mass1 = planet1.mass;
    vec pos2 = planet2.position;
    double mass2 = planet2.mass;

    double r = norm(pos2-pos1); // Relative distance between the objects.

    // The gravitational force, now with a different power (beta) in the
    // inverse law:
    double forceStrength = (G*mass1*mass2)/pow(r,beta);
    vec forceDirection = (pos2-pos1)/r;
    vec forceVector = forceStrength * forceDirection;
    return forceVector;
}

vec gForceGenRelCorr(Planet planet1, Planet planet2, double G){
    // This returns the gravitational force *on* planet 1 *from* planet 2. Object 1
    // is pulled towards object 2.
    double c = kmPerSec_to_auPerYear(299792.458);
    vec pos1 = planet1.position;
    double mass1 = planet1.mass;
    vec pos2 = planet2.position;
    double mass2 = planet2.mass;
    double angMom1 = planet1.angularMomentum()/planet1.mass; 

    double r = norm(pos2 - pos1); // Relative distance between the objects.

    double forceStrength = (1 + 3*angMom1*(1/(r*r*c*c)))*(G*mass1*mass2)*(1/(r*r));   // Newton's gravitational law.

    vec forceDirection = (pos2-pos1)/norm(pos2-pos1);   // This vector points *from*
    // object 1 and *towards* object 2, meaning that object 1 is influenced by object 2.
    vec forceVector = forceStrength * forceDirection;
    //cout << forceVector.t() << endl;
    return forceVector;
}

double potentialEnergy(Planet current, Planet other, double G){
    return -G*other.mass*current.mass/norm(current.position - other.position);
}

void writeMatrixToFile(mat results, string fileName, string directory){
    // Write the results (an Nx7 matrix) from an ODE solver to
    // a text file with 7 columns.
    // fileName: The full name of the file, e.g. "data.txt".
    // directory: Specify the directory where the file is to be saved. E.g.
    // "../results/3a_earth_sun_system/" (include the final slash).
    
    ofstream ofile;
    string filePath = directory + fileName;

    // Save matrix in CSV format with a header:
    field<string> header(results.n_cols);
    header(0) = "t"; header(1) = "x"; header(2) = "y"; header(3) = "z";
    header(4) = "vx"; header(5) = "vy"; header(6) = "vz";
    //results.save(csv_name("results.csv", header));
    results.save(csv_name(filePath, header));
}

void writeGeneralMatrixToCSV(mat results, field<string> columnLabels, string filename, string directory){
    // columnLabels contains the labels of the columns, e.g. "t", "x", "y", "z" or "t", "L".
    // It should have the same amount of elements as 'results' has columns.
    ofstream ofile;
    string filePath = directory + filename;

    // Save matrix in CSV format with the column labels in the header:
    //results.save(csv_name("results.csv", header));
    results.save(csv_name(filePath, columnLabels));
}

void writeSolarSystemToFiles(mat resultsAllPlanets, int nTimesteps, int nPlanets, field<string> planetNames, string directory){
    // This function writes the results matrix from Solver::run_velocityVerletForceType()
    // for multiple planets to individual .csv files, one for each planet.
    // The results matrix is constructed such that the first N rows (after the header
    // in the first row) contains the positions and velocities for the first planet (the
    // Sun not included), the next N rows contains the positions and velocities for the
    // second planet, and so on.
    // planetNames: A vector containing the names of the different planets. They should
    // be in the same order as the planets are added to the Solver object (the solar
    // system). This is going to be included in their respective file names.
    // nTimesteps: The number of timesteps in the simulation.
    // nPlanets: The number of planets.

    // Create one textfile for each planet:
    for (int i=0; i<=nPlanets; i++){
        // File name:
        string planetName = planetNames[i];
        string fileName = planetName + ".csv";

        int startIdx = 0 + i*nTimesteps; // Index of the first timestep
        // of the current planet in the resultsAllPlanets matrix.
        // Get the sub-matrix that corresponds to the results of the current
        // planet (nTimesteps rows and all 7 columns):
        mat resultsThisPlanet = resultsAllPlanets(span(startIdx, startIdx+(nTimesteps-1)), span(0, 6));
        //resultsThisPlanet.print(planetName);

        // Write the sub-matrix for the current planet to file:
        writeMatrixToFile(resultsThisPlanet, fileName, directory);
    }
}

double kmPerSec_to_auPerYear(double speed_kmPerSec){
    double oneKmPerSec_in_auPerYear = 0.21094502111897098; // This is 1 km/sec in units au/yr.
    return speed_kmPerSec * oneKmPerSec_in_auPerYear;
}

double auPerYear_to_kmPerSec(double speed_auPerYear){
    double oneAuPerYear_in_kmPerSec = 4.74057171245587; // This is 1 au/yr in units km/sec
    return speed_auPerYear * oneAuPerYear_in_kmPerSec;
}

double get_earth_mass(){
    //double m_S_SI = 1.989e30; // kg
    //double m_E_SI = 5.972e24; // kg
    //double m_E = m_E_SI/m_S_SI;
    double m_E = 3.002463426e-6; // Earth's mass in solar masses
    return m_E;
}

double get_mercury_mass(){
    return 3.285e23/1.989e30;
}

vec initial_earth_velocity(vec initialPosition){
    // Return the earth initial velocity from position. 
    
    vec omegaDirection = vec("0 0 1");
    double R_E = 1.;    // Earth orbit radius (1 AU).
    double m_S = 1.; // Solar mass in units of solar masses.
    double T_E = 1.; // Orbit period of Earth: 1 year.
    double omega_E = (2*M_PI)/T_E; // Angular frequency of Earth's orbit (radians).
    double v_E = (2*M_PI*R_E)/T_E;  // Orbital speed of Earth (approximately constant along the
    // entire orbit since the orbit is approximately circular).

    vec v_E_dir = cross(omegaDirection, initialPosition) / norm(cross(omegaDirection, initialPosition)); // The direction of the orbital velocity
    // is perpendicular to both the angular momentum and the position in the orbit.
    vec initialVelocity = v_E * v_E_dir;
    return initialVelocity;
}

mat run_forwardEuler(double tFinal, double dt, double G){
    // This function returns an Nx7 matrix where the columns contain the solution
    // data points (t, x, y, z, vx, vy, vz) (positions and velocities and their corresponding
    // times) gotten from forward Euler. N is the number of timesteps.

    // Here: Choose the x-y plane as the orbit plane of Earth. Then the
    // angular momentum vector (omega) of Earth in orbit will point in the positive
    // z direction.
    vec omegaDirection = vec("0 0 1");

    int N = round(tFinal/dt);   // Number of timesteps.
    cout << "Number of timesteps: N = " << N << endl;
    vec tList = vec(N);
    for(int i=0; i<=N-1; i++){tList(i) = i*dt;}

    // Initial position of Earth:
    vec initialPosition = vec("1 0 0");

    vec sunPosition = vec("0 0 0"); // The sun is at the center and approximately
    // doesn't move because of its relatively large mass.

    // Masses of Sun and Earth in SI units (kg):
    double m_S_SI = 1.989e30;
    double m_E_SI = 5.972e24;
    double R_E = 1.;    // Earth orbit radius (1 AU).
    
    double m_S = 1.; // Solar mass in units of solar masses.
    //double m_S = m_S_SI;
    double m_E = m_E_SI/m_S_SI;     // Earth mass in units of solar masses.
    //double m_E = m_E_SI;

    double T_E = 1.; // Orbit period of Earth: 1 year.
    double omega_E = (2*M_PI)/T_E; // Angular frequency of Earth's orbit (radians).
    double v_E = (2*M_PI*R_E)/T_E;  // Orbital speed of Earth (approximately constant along the
    // entire orbit since the orbit is approximately circular).
    vec v_E_dir = cross(omegaDirection, initialPosition) / norm(cross(omegaDirection, initialPosition)); // The direction of the orbital velocity
    // is perpendicular to both the angular momentum and the position in the orbit.
    vec initialVelocity = v_E * v_E_dir;

    // Initialize positions and velocities:
    vec position = initialPosition; // <-- An example of the starting position of Earth;
    // 1 AU away from the Sun in the x direction.
    vec velocity = initialVelocity;
    
    
    // Initialize acceleration vector (from Newton's 2nd law):
    vec FVec = gForceVector(G, m_E, m_S, initialPosition, sunPosition);
    vec initialAcc = FVec/m_E; // The acceleration of Earth (from Newton's 2nd law).
    vec acceleration = initialAcc;

    mat results = mat(N, 7); // Columns: t, x, y, z, vx, vy, vz
    results.col(0) = tList; // Insert the time list in the results matrix.
    results(0, span(1,3)) = initialPosition.t();
    results(0, span(4,6)) = initialVelocity.t();

    // Perform the forward Euler algorithm:
    for(int i=1; i<=N-1; i++){
        // Update the position using the previous position and velocity:
        position = position + dt*velocity;

        // Update the velocity using the acceleration:
        velocity = velocity + dt*acceleration;

        results(i, span(1,3)) = position.t(); // Transpose the vectors in order
        // to paste them into the row.
        results(i, span(4,6)) = velocity.t();

        // Update the acceleration for the next time step (to update
        // the velocity):
        
        FVec = gForceVector(G, m_E, m_S, position, sunPosition);
        //FVec = gForceVector(G, m_E, m_S, position, sunPosition) * (1/2.976e-19); // Debugging (the factor is to check if G must be of other units.)
        acceleration = FVec/m_E;
    }

    return results;
}

mat forwardEuler(double tFinal, double dt, double m_SI, vec initialPosition, vec initialVelocity, double G){
    // This is the same as run_forwardEuler(), but more flexible since it has more input arguments.
    // This function returns an Nx7 matrix where the columns contain the solution
    // data points (t, x, y, z, vx, vy, vz) (positions and velocities and their corresponding
    // times) gotten from forward Euler. N is the number of timestep
    // initialPosition: 3D vector, the initial position of the planet in astronomical units.
    // initialVelocity: 3D vector, the initial velocity of the planet in au/yr.
    // m_SI: Mass of the planet in units kg. Earth mass: 5.972e24 kg.

    int N = round(tFinal/dt);   // Number of timesteps.
    cout << "Number of timesteps: N = " << N << endl;
    vec tList = vec(N);
    for(int i=0; i<=N-1; i++){tList(i) = i*dt;}

    vec sunPosition = vec("0 0 0"); // The sun is at the center and approximately
    // doesn't move because of its relatively large mass.

    // Masses of Sun and Earth in SI units (kg):
    double m_S_SI = 1.989e30;
    
    double m_S = 1.; // Solar mass in units of solar masses.
    double m = m_SI/m_S_SI;     // Planet mass in units of solar masses.

    // Initialize positions and velocities:
    vec position = initialPosition; // <-- An example of the starting position of Earth;
    // 1 AU away from the Sun in the x direction.
    vec velocity = initialVelocity;
    
    // Initialize acceleration vector (from Newton's 2nd law):
    vec FVec = gForceVector(G, m, m_S, initialPosition, sunPosition);
    vec initialAcc = FVec/m; // The acceleration of Earth (from Newton's 2nd law).
    vec acceleration = initialAcc;

    mat results = mat(N, 7); // Columns: t, x, y, z, vx, vy, vz
    results.col(0) = tList; // Insert the time list in the results matrix.
    results(0, span(1,3)) = initialPosition.t();
    results(0, span(4,6)) = initialVelocity.t();

    // Perform the forward Euler algorithm:
    for(int i=1; i<=N-1; i++){
        // Update the position using the previous position and velocity:
        position = position + dt*velocity;

        // Update the velocity using the acceleration:
        velocity = velocity + dt*acceleration;

        results(i, span(1,3)) = position.t(); // Transpose the vectors in order
        // to paste them into the row.
        results(i, span(4,6)) = velocity.t();

        // Update the acceleration for the next time step (to update
        // the velocity): 
        FVec = gForceVector(G, m, m_S, position, sunPosition);
        //FVec = gForceVector(G, m_E, m_S, position, sunPosition) * (1/2.976e-19); // Debugging (the factor is to check if G must be of other units.)
        acceleration = FVec/m;
    }

    return results;
}

mat run_velocityVerlet(double tFinal, double dt, double G){
    vec omegaDirection = vec("0 0 1");

    int N = round(tFinal/dt);   // Number of timesteps.
    cout << "Number of timesteps: N = " << N << endl;
    vec tList = vec(N);
    for(int i=0; i<=N-1; i++){tList(i) = i*dt;}

    // Initial position of Earth:
    vec initialPosition = vec("1 0 0");

    vec sunPosition = vec("0 0 0"); // The sun is at the center and approximately
    // doesn't move because of its relatively large mass.

    // Masses of Sun and Earth in SI units (kg):
    double m_S_SI = 1.989e30;
    double m_E_SI = 5.972e24;
    double R_E = 1.;    // Earth orbit radius (1 AU).
    
    double m_S = 1.; // Solar mass in units of solar masses.
    double m_E = m_E_SI/m_S_SI;     // Earth mass in units of solar masses.

    double T_E = 1.; // Orbit period of Earth: 1 year.
    double omega_E = (2*M_PI)/T_E; // Angular frequency of Earth's orbit (radians).
    double v_E = (2*M_PI*R_E)/T_E;  // Orbital speed of Earth (approximately constant along the
    // entire orbit since the orbit is approximately circular).
    vec v_E_dir = cross(omegaDirection, initialPosition) / norm(cross(omegaDirection, initialPosition)); // The direction of the orbital velocity
    // is perpendicular to both the angular momentum and the position in the orbit.
    vec initialVelocity = v_E * v_E_dir;

    // Initialize positions and velocities:
    vec position = initialPosition; // <-- An example of the starting position of Earth;
    // 1 AU away from the Sun in the x direction.
    vec velocity = initialVelocity;

    // Initialize acceleration vector (from Newton's 2nd law):
    vec FVec = gForceVector(G, m_E, m_S, initialPosition, sunPosition);
    vec initialAcc = FVec/m_E; // The acceleration of Earth (from Newton's 2nd law).
    vec acceleration = initialAcc;

    mat results = mat(N, 7); // Columns: t, x, y, z, vx, vy, vz
    results.col(0) = tList; // Insert the time list in the results matrix.
    results(0, span(1,3)) = initialPosition.t();
    results(0, span(4,6)) = initialVelocity.t();

    vec previous_acceleration;
    // Perform the velocity Verlet algorithm:
    for(int i=1; i<=N-1; i++){
        previous_acceleration = acceleration;

        // Evaluate the current position, acceleration and velocity.
        position += dt*velocity + 0.5*dt*dt*previous_acceleration;
        acceleration = gForceVector(G, m_E, m_S, position, sunPosition)/m_E;
        velocity += 0.5*dt*(acceleration + previous_acceleration);

        results(i, span(1,3)) = position.t(); // Transpose the vectors in order
        results(i, span(4,6)) = velocity.t();

    }
    return results;
}

mat velocityVerlet(double tFinal, double dt, double m_SI, vec initialPosition, vec initialVelocity, double G){
    int N = round(tFinal/dt);   // Number of timesteps.
    cout << "Number of timesteps: N = " << N << endl;
    vec tList = vec(N);
    for(int i=0; i<=N-1; i++){tList(i) = i*dt;}

    vec sunPosition = vec("0 0 0"); // The sun is at the center and approximately
    // doesn't move because of its relatively large mass.

    // Masses of Sun and Earth in SI units (kg):
    double m_S_SI = 1.989e30;
    
    double m_S = 1.; // Solar mass in units of solar masses.
    double m = m_SI/m_S_SI;     // Planet mass in units of solar masses.

    // Initialize positions and velocities:
    vec position = initialPosition; // <-- An example of the starting position of Earth;
    // 1 AU away from the Sun in the x direction.
    vec velocity = initialVelocity;
    
    // Initialize acceleration vector (from Newton's 2nd law):
    vec FVec = gForceVector(G, m, m_S, initialPosition, sunPosition);
    vec initialAcc = FVec/m; // The acceleration of Earth (from Newton's 2nd law).
    vec acceleration = initialAcc;

    mat results = mat(N, 7); // Columns: t, x, y, z, vx, vy, vz
    results.col(0) = tList; // Insert the time list in the results matrix.
    results(0, span(1,3)) = initialPosition.t();
    results(0, span(4,6)) = initialVelocity.t();

    vec previous_acceleration;
    // Perform the velocity Verlet algorithm:
    for(int i=1; i<=N-1; i++){
        previous_acceleration = acceleration;

        // Evaluate the current position, acceleration and velocity.
        position += dt*velocity + 0.5*dt*dt*previous_acceleration;
        acceleration = gForceVector(G, m, m_S, position, sunPosition)/m;
        velocity += 0.5*dt*(acceleration + previous_acceleration);

        results(i, span(1,3)) = position.t(); // Transpose the vectors in order
        results(i, span(4,6)) = velocity.t();

    }
    return results;
}

void task_3a_forwardEuler(double G){
    // This runs problem 3a with the forward Euler algorithm.
    cout << "Running forward Euler:\n";
    double tFinal = 10;
    double dt = 1e-4;
    cout << "tFinal: " << tFinal << " (years) \ndt: " << dt << endl;

    mat resultsEuler = run_forwardEuler(tFinal, dt, G);
    string filename = "earth_sun_euler.csv";
    string directory = "../results/3a_earth_sun_system/";
    writeMatrixToFile(resultsEuler, filename, directory);
}

void task_3a_velocityVerlet(double G){
    // Runs the problem 3a with velocity verlet algorithm.
    cout << "Running velocity Verlet:\n";
    double tFinal = 10;
    double dt = 1e-4;
    cout << "tFinal: " << tFinal << " (years) \ndt: " << dt << endl;

    // Get the results:
    mat resultsVerlet = run_velocityVerlet(tFinal, dt, G);

    string filename = "earth_sun_verlet.csv";
    string directory = "../results/3a_earth_sun_system/";
    writeMatrixToFile(resultsVerlet, filename, directory); 
}

void task_3b_velocityVerlet(double G){
    // Runs object oriented velocity Verlet. 

    double dt = 1e-3;
    //double tFinal = 10;
    int N = 1e4;
    double tFinal = dt*N;

    // Initial position and velocity of Earth.
    vec initialPosition = vec("1 0 0");
    vec initialVelocity = initial_earth_velocity(initialPosition);

    // Initial pos and vel of Sun.
    vec sunPosition = vec("0 0 0"); 
    vec sunVelocity = vec("0 0 0");
    
    double m_S = 1.0;
    Planet sun;
    sun.init(m_S, sunPosition, sunVelocity);    

    double m_E = get_earth_mass();
    Planet earth;
    earth.init(m_E, initialPosition, initialVelocity); 

    Solver my_solver;
    my_solver.init(N);
    my_solver.add(sun);
    my_solver.add(earth);

    mat resultsVerlet = my_solver.run_velocityVerletForceType(0, tFinal, dt, G);
    string filename = "earth_sun_verlet.csv";
    string directory = "../results/3b_earth_sun_system/";
    writeMatrixToFile(resultsVerlet, filename, directory); 

    mat momEnergyMatrix = my_solver.get_angMomentum_energy_mat();
    string filename1 = "earth_sun_energy.csv";
    momEnergyMatrix.save(csv_name(directory + filename1));
}

void task_3e_force(double G){
    // Runs object oriented velocity Verlet with varying versions of
    // the gravitational force (determined by beta).

    // The directory where the results will be saved:
    string directory = "../results/3e_force/";

    // Choose which values of beta to run through:
    /*
    double betaMin = 2; double betaMax = 3; // Limit values of beta
    int betaListLength = 6; // Number of elements in the beta list
    double betaStepSize = (betaMax-betaMin)/(betaListLength-1); // Step size in the beta list
    cout << "betaStepSize: " << betaStepSize << endl;
    vec betaList = vec(betaListLength);
    betaList(0) = betaMin;
    betaList(betaListLength-1) = betaMax;
    // Fill betaList:
    for (int i=1; i<=betaListLength-2; i++){betaList(i) = betaList(0)+i*betaStepSize;}
    */
    /*
    double betaMin = 2.0; // Limit values of beta
    int betaListLength = 6; // Number of elements in the beta list
    double betaStepSize = 0.2; // Step size in the beta list
    vec betaList = vec(betaListLength);
    betaList(0) = betaMin;
    // Fill betaList:
    for (int i=1; i<=betaListLength-1; i++){betaList(i) = betaList(0)+i*betaStepSize;}
    */
    //vec betaList = vec("2.0 2.2 2.4 2.6 2.8 2.9 3.0");
    vec betaList = vec("2.0 2.2 2.4 2.6 2.8 2.9 3.0");
    double betaListLength = betaList.n_elem;
    //cout << "betaListLength" << betaListLength << endl;

    // Save the beta values in a file so that the python script
    // can get the beta values:
    betaList.save(csv_name(directory + "betaList.csv"));

    // Timestep and iterations:
    double dt = 1e-4;
    double tFinal = 5; int N = round(tFinal/dt);
    //int N = 500; double tFinal = dt*N;

    double m_S = 1.0;
    double m_E = get_earth_mass();

    // Initial position and velocity of Earth:
    vec initialPosition_E = vec("1 0 0");
    // If doing the first part of 3e, choose circular orbit velocity:
    //vec initialVelocity_E = initial_earth_velocity(initialPosition_E);
    // If doing the second part of 3e, choose speed 5 AU/yr:
    vec initialVelocity_E = vec("0 5 0");

    // Position and velocity of the Sun:
    vec sunPosition = vec("0 0 0"); 
    vec sunVelocity = vec("0 0 0");

    // For each beta value, get and print the solar system evolution:
    //for (auto beta : betaList){
    for (int i=0; i<=betaListLength-1; i++){
        // The planets and the solver must be re-initialized for each value of beta
        // since we're starting over again for each beta value.
        double beta = betaList(i);

        Planet sun;
        sun.init(m_S, sunPosition, sunVelocity);
        Planet earth;
        earth.init(m_E, initialPosition_E, initialVelocity_E); 

        Solver my_solver;
        my_solver.init(N);
        my_solver.add(sun);
        my_solver.add(earth);

        //cout << "beta: " << beta << endl;
        mat results = my_solver.run_velocityVerletBeta(tFinal, dt, beta, G);

        // Save one results file for each beta value:
        stringstream betaStringstream;  // The numerical value of beta as a string
        betaStringstream << setprecision(2) << fixed << beta;
        string betaStr = betaStringstream.str();
        // File name:
        string filename = "3e_force_beta" + betaStr + ".csv";
        //cout << "filename (beta): " << filename << endl;
        // Write the results matrix to the results file:
        writeMatrixToFile(results, filename, directory);
    }
}

void task_3f_escape_velocity(double initialSpeed_kmPerSec, double G){
    // This does essentially the same as task_3a_velocityVerlet(), but
    // customized for task 3f.
    // initialSpeed: Initial orbit speed in km/s.

    double initialSpeed_auPerYear = kmPerSec_to_auPerYear(initialSpeed_kmPerSec);
    cout << "Initial orbit speed in km/sec: " << initialSpeed_kmPerSec << " km/sec = " 
    << initialSpeed_auPerYear << " au/yr.\n";

    double tFinal = 500;
    double dt = 1e-3;

    vec initialPosition = vec("1 0 0");
    vec initialVelocityDirection = vec("0 1 0");
    vec initialVelocity = initialSpeed_auPerYear * initialVelocityDirection;

    double m_SI = 5.972e24; // The mass of Earth (kg) (but the mass doesn't matter,
    // no pun intended.)

    cout << "Running velocity Verlet (testing escape velocity):\n";
    cout << "tFinal: " << tFinal << " (years) \ndt: " << dt << endl;

    mat resultsVerlet = velocityVerlet(tFinal, dt, m_SI, initialPosition, initialVelocity, G);

    string filename = "test_escape_verlet.csv";
    string directory = "../results/3f_escape_velocity/";
    writeMatrixToFile(resultsVerlet, filename, directory);
}

void task_3g_three_body(double G){
    // Runs object oriented velocity Verlet. 
    double dt = 1e-4;
    //double tFinal = 5; int N = round(tFinal/dt);
    int N = 50; double tFinal = dt*N;
    
    
    // Initial positions and velocities:
    // Sun:
    vec sunPosition = vec("0 0 0");
    vec sunVelocity = vec("0 0 0");
    // Earth:
    vec initialPosition = vec("1 0 0");
    vec initialVelocity = initial_earth_velocity(initialPosition);
    // Jupiter: (Retrieved from NASA webpage)
    vec initPosJupiter = vec("5.2044 0 0"); // au
    vec initVelJupiter = vec("0 2.75522 0"); // au/yr

    // Initialize planets (the Sun must be initialized first):    
    Planet sun;
    double m_S = 1.0;
    sun.init(m_S, sunPosition, sunVelocity);    

    Planet earth;
    double m_E = get_earth_mass();
    earth.init(m_E, initialPosition, initialVelocity);

    Planet jupiter;
    double massMultiplier = 1; // How much we want to multiply Jupiter's mass by
    double m_J = 9.545536837e-4 * massMultiplier; // in solar masses
    jupiter.init(m_J, initPosJupiter, initVelJupiter);

    Solver my_solver;
    my_solver.init(N);
    my_solver.add(sun);
    my_solver.add(earth);
    my_solver.add(jupiter);

    // Get the number of planets (excluding the Sun):
    int nPlanets = my_solver.get_total_planets()-1;

    // The results matrix contains the results for all planets:
    mat resultsAllPlanets = my_solver.run_velocityVerletForceType(0, tFinal, dt, G);

    // The results directory:
    string directory = "../results/3g_three_body/";

    // Store the planet names in a field:
    field<string> planetNames(nPlanets);
    planetNames(0) = "earth"; planetNames(1) = "jupiter";
    // Save the planet names to a file so that the python script can get them:
    planetNames.save(directory + "planet_names.csv");

    // Write the results for the planets to files:
    writeSolarSystemToFiles(resultsAllPlanets, N, nPlanets, planetNames, directory);
}

void task_3h_solar_system(double G){
    // Runs object oriented velocity Verlet. 
    double dt = 1e-4;
    double tFinal = 5; int N = round(tFinal/dt);
    //int N = 1000; double tFinal = dt*N;
    
    // Planets: Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
    
    // Initial positions and velocities:
    // Sun:
    vec sunPosition = vec("0 0 0");
    vec sunVelocity = vec("0 0 0");
    // Mercury:
    vec initPosMercury = vec("0.34439581 -0.15170017 -0.04498414");
    vec initVelMercury = vec("2.21697462 9.83032311 0.59991024");
    // Venus:
    vec initPosVenus = vec("-0.18184097  0.70330564  0.01979342");
    vec initVelVenus = vec("-7.19107188 -1.84968568  0.38953739");
    // Earth:
    vec initPosEarth = vec("9.26247918e-01 3.60984141e-01 6.66364627e-05");
    vec initVelEarth = vec("-2.33799168e+00  5.85066441e+00 -2.18690746e-04");
    // Mars:
    vec initPosMars = vec("1.31769404  0.5095586  -0.02184232");
    vec initVelMars = vec("-1.62286702  5.21232688  0.14909549");
    // Jupiter:
    vec initPosJupiter = vec("2.54375857 -4.4368376  -0.0385057"); // au
    vec initVelJupiter = vec("2.35721704  1.50155279 -0.05896163"); // au/yr
    // Saturn:
    vec initPosSaturn = vec("5.13508534 -8.56941586 -0.05543229");
    vec initVelSaturn = vec("1.63444679  1.04237547 -0.08332824");
    // Uranus:
    vec initPosUranus = vec("15.53976957 12.23945879 -0.15586173");
    vec initVelUranus = vec("-0.89936374  1.0616034   0.01560822");
    // Neptune:
    vec initPosNeptune = vec("29.41117938 -5.47128787 -0.56513995");
    vec initVelNeptune = vec("0.20218953  1.13396172 -0.02815956");
    // Pluto:
    vec initPosPluto = vec("13.82177646 -31.20040462  -0.65945608");
    vec initVelPluto = vec("1.07229549  0.21866022 -0.33001162");

    // Initialize the planets (the Sun must be initialized first):    
    Planet sun;
    double m_S = 1.0;
    sun.init(m_S, sunPosition, sunVelocity);  
    
    Planet mercury;
    double m_Mercury = 1.6600955494091023e-07;
    mercury.init(m_Mercury, initPosMercury, initVelMercury);
    
    Planet venus;
    double m_V = 2.447824993713855e-06;
    venus.init(m_V, initPosVenus, initVelVenus);
    
    Planet earth;
    double m_E = get_earth_mass();
    earth.init(m_E, initPosEarth, initVelEarth);
    
    Planet mars;
    double m_Mars = 3.2271058586874533e-07;
    mars.init(m_Mars, initPosMars, initVelMars);
    
    Planet jupiter;
    double m_J = 9.545536837e-4; // in solar masses
    jupiter.init(m_J, initPosJupiter, initVelJupiter);
    
    Planet saturn;
    double m_Saturn = 0.00028581342720643705;
    saturn.init(m_Saturn, initPosSaturn, initVelSaturn);

    Planet uranus;
    double m_U = 4.365602212723158e-05;
    uranus.init(m_U, initPosUranus, initVelUranus);

    Planet neptune;
    double m_N = 5.150264018104099e-05;
    neptune.init(m_N, initPosNeptune, initVelNeptune);

    Planet pluto;
    double m_P = 6.552677897913e-09;
    pluto.init(m_P, initPosPluto, initVelPluto);

    Solver my_solver;
    my_solver.init(N);
    my_solver.add(sun);
    my_solver.add(mercury);
    my_solver.add(venus);
    my_solver.add(earth);
    my_solver.add(mars);
    my_solver.add(jupiter);
    my_solver.add(saturn);
    my_solver.add(uranus);
    my_solver.add(neptune);
    my_solver.add(pluto);

    // Get the number of planets (excluding the Sun):
    int nPlanets = my_solver.get_total_planets()-1;
    cout << "nPlanets: " << nPlanets << endl;

    // The results matrix contains the results for all planets:
    mat resultsAllPlanets = my_solver.run_velocityVerletForceType(0, tFinal, dt, G);
    //mat resultsAllPlanets = my_solver.run_velocityVerletBeta(tFinal, dt, 2.0, G);

    // The results directory:
    string directory = "../results/3h_solar_system/";

    // Store the planet names in a field:
    field<string> planetNames(nPlanets);
    planetNames(0) = "mercury";
    planetNames(1) = "venus";
    planetNames(2) = "earth";
    planetNames(3) = "mars";
    planetNames(4) = "jupiter";
    planetNames(5) = "saturn";
    planetNames(6) = "uranus";
    planetNames(7) = "neptune";
    planetNames(8) = "pluto";

    // Save the planet names to a file so that the python script can get them:
    planetNames.save(directory + "planet_names.csv");

    // Write the results for the planets to files:
    writeSolarSystemToFiles(resultsAllPlanets, N, nPlanets, planetNames, directory);
    cout << "debugging: end of utils.cpp\n";
}


void task_3i_mercury_precession(double G){
    // Runs object oriented velocity Verlet. 

    double dt = 1e-5;
    double tFinal = 100;
    int N = round(tFinal/dt);
    cout << N << endl;

    // Initial position and velocity of Mercury.
    vec initialPosition = vec("0.307499 0 0");
    vec initialVelocity = vec("0 12.44 0");

    // Initial pos and vel of Sun.
    vec sunPosition = vec("0 0 0"); 
    vec sunVelocity = vec("0 0 0");
    
    double m_S = 1.0;
    Planet sun;
    sun.init(m_S, sunPosition, sunVelocity);    
    
    double m_E = get_mercury_mass();
    Planet mercury;
    mercury.init(m_E, initialPosition, initialVelocity); 

    Solver my_solver;
    my_solver.init(N);
    my_solver.add(sun);
    my_solver.add(mercury);

    // 3D matrix instead? One layer for each matrix? (from run_velocityVerlet)
    mat resultsVerlet = my_solver.run_velocityVerletForceType(1, tFinal, dt, G);
    string filename = "mercury_sun_verlet.csv";
    string directory = "../results/3i_mercury_precession/";
    writeMatrixToFile(resultsVerlet, filename, directory); 

    mat momEnergyMatrix = my_solver.get_angMomentum_energy_mat();
    string filename1 = "mercury_sun_energy.csv";
    momEnergyMatrix.save(csv_name(directory + filename1));

    vector<Planet> all_planets = my_solver.get_all_planets();
    mat peri = my_solver.perihelion_mat_solver;
    string filename_peri = "perihelion.csv";

    field<string> header(peri.n_cols);
    header(0) = "x"; header(1) = "y"; header(2) = "z";
    writeGeneralMatrixToCSV(peri, header, filename_peri, directory);
    
}