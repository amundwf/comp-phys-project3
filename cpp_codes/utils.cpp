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

/*
void run_forwardEuler(){
    // Here: Choose the x-y plane as the orbit plane of Earth. Then the
    // angular momentum vector (omega) of Earth in orbit will point in the positive
    // z direction.
    vec omegaDirection = vec("0 0 1");

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
    
    double G = 6.67e-11; // N m^2/kg^2. Gravitational constant.
    // Initialize acceleration (from Newton's 2nd law):
    double Fx, Fy, Fz;
    vec F = gForceVector(G, m_E, m_S, initialPosition, sunPosition);

    // Initialize forces
    double Fx,Fy,Fz,Fxnew,Fynew,Fznew; // Forces in each dimension
}
*/

/*
void run_velocityVerlet(){

}
*/