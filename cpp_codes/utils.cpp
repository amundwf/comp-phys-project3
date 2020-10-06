#include "utils.hpp" 
#include <iostream>
#include <cmath>
#include <armadillo>

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
    return forceVector;
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

    //double GM_s = 4*M_PI*M_PI; // AU^3/yr^2 (from centripetal acceleration)

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

void task_3a_forwardEuler(double G){
    // This runs problem 3a with the forward Euler algorithm.
    double tFinal = 10;
    double dt = 1e-4;
    mat resultsEuler = run_forwardEuler(tFinal, dt, G);
    string filename = "earth_sun_euler.csv";
    string directory = "../results/3a_earth_sun_system/";
    writeMatrixToFile(resultsEuler, filename, directory);
}


void task_3a_velocityVerlet(double G){
    // Runs the problem 3a with velocity verlet algorithm.
    double tFinal = 10;
    double dt = 1e-4;
    mat resultsVerlet = run_velocityVerlet(tFinal, dt, G);
    string filename = "earth_sun_verlet.csv";
    string directory = "../results/3a_earth_sun_system/";
    writeMatrixToFile(resultsVerlet, filename, directory); 
}