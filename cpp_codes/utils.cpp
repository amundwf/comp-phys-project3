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
    // NB: If r=1, then pow(r,beta) = r^beta = 1^beta = 1 for all values of
    // beta!

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

void writeMatrixToFile(mat results, string filename, string directory){
    // Write the results (an Nx7 matrix) from an ODE solver to
    // a text file with 7 columns.
    // filename: The full name of the file, e.g. "data.txt".
    // directory: Specify the directory where the file is to be saved. E.g.
    // "../results/3a_earth_sun_system/" (include the final slash).
    
    ofstream ofile;
    string filePath = directory + filename;

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

    // 3D matrix instead? One layer for each matrix? (from run_velocityVerlet)
    mat resultsVerlet = my_solver.run_velocityVerlet(tFinal, dt, G);
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
    double betaMin = 2.0; // Limit values of beta
    int betaListLength = 6; // Number of elements in the beta list
    double betaStepSize = 0.2; // Step size in the beta list
    vec betaList = vec(betaListLength);
    betaList(0) = betaMin;
    // Fill betaList:
    for (int i=1; i<=betaListLength-1; i++){betaList(i) = betaList(0)+i*betaStepSize;}
    //vec betaList = vec("1.0, 2.0, 2.5, 3.0, 4.0");

    // Save the beta values in a file so that the python script
    // can get the beta values:
    betaList.save(csv_name(directory + "betaList.csv"));

    // Timestep and iterations:
    double dt = 1e-4;
    double tFinal = 3; int N = round(tFinal/dt);
    //int N = 500; double tFinal = dt*N;

    double m_S = 1.0;
    double m_E = get_earth_mass();

    // Initial position and velocity of Earth:
    vec initialPosition_E = vec("1 0 0");
    // If doing the first part of 3e, choose circular orbit velocity:
    //vec initialVelocity_E = initial_earth_velocity(initialPosition_E);
    // If doing the second part of 3e, choose speed 5 AU/yr:
    vec initialVelocity_E = vec("0 5 0");

    // initialVelocity_E = initial_earth_velocity_beta(initialPosition);
    // ^ Does the velocity need to be different? Do we want to keep the same
    // velocity or do we want to change the velocity to match the new centripetal
    // acceleration from the new force law?
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
        betaStringstream << setprecision(1) << fixed << beta;
        string betaStr = betaStringstream.str();
        // File name:
        string filename = "3e_force_beta" + betaStr + ".csv";
        //cout << "filename (beta): " << filename << endl;
        // Write the results matrix to the results file:
        writeMatrixToFile(results, filename, directory);
    }

    // Write the number of beta values to a file so that the python
    // script for this task can automatically read the number of
    // beta values?:
    // writeToFile('numBetaValues', betaList.n_elem); // <-- or something like that
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
    double tFinal = 0.75;
    //int N = 5000; double tFinal = dt*N;
    int N = round(tFinal/dt);
    
    // Initial position and velocity of Earth.
    vec initialPosition = vec("1 0 0");
    vec initialVelocity = initial_earth_velocity(initialPosition);

    // Initial pos and vel of Sun.
    vec sunPosition = vec("0 0 0"); 
    vec sunVelocity = vec("0 0 0");
    
    initialVelocity.t().print("initialVelocity: ");
    
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

    // 3D matrix instead? One layer for each matrix? (from run_velocityVerlet)
    mat resultsVerlet = my_solver.run_velocityVerlet(tFinal, dt, G);
    // my_solver.run_velocityVerlet(tFinal, dt, G); ? 
    //resultsVerlet.print("resultsVerlet");

    string filename = "three_body_verlet.csv";
    string directory = "../results/3g_three_body/";
    writeMatrixToFile(resultsVerlet, filename, directory); 
}

void task_3i_mercury_precession(double G){
    // Runs object oriented velocity Verlet. 

    double dt = 1e-3;
    double tFinal = 100;
    int N = round(tFinal/dt);
    cout << N << endl;

    // Initial position and velocity of Mercury.
    vec initialPosition = vec("1 0 0");
    vec initialVelocity = initial_earth_velocity(initialPosition);
    //vec initialPosition = vec("0.307499 0 0");
    //vec initialVelocity = vec("12.44 0 0");

    // Initial pos and vel of Sun.
    vec sunPosition = vec("0 0 0"); 
    vec sunVelocity = vec("0 0 0");
    
    double m_S = 1.0;
    Planet sun;
    sun.init(m_S, sunPosition, sunVelocity);    
    /*
    double m_E = get_mercury_mass();
    Planet mercury;
    mercury.init(m_E, initialPosition, initialVelocity); 
    */
    double m_E = get_earth_mass();
    Planet earth;
    earth.init(m_E, initialPosition, initialVelocity); 

    Solver my_solver;
    my_solver.init(N);
    my_solver.add(sun);
    my_solver.add(earth);
    //my_solver.add(mercury);

    // 3D matrix instead? One layer for each matrix? (from run_velocityVerlet)
    mat resultsVerlet = my_solver.run_velocityVerletForceType(1, tFinal, dt, G);
    string filename = "mercury_sun_verlet.csv";
    string directory = "../results/3i_mercury_precession/";
    writeMatrixToFile(resultsVerlet, filename, directory); 

    mat momEnergyMatrix = my_solver.get_angMomentum_energy_mat();
    string filename1 = "mercury_sun_energy.csv";
    momEnergyMatrix.save(csv_name(directory + filename1));
}