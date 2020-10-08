#include "planet.hpp"
#include <armadillo>

using namespace arma;
using namespace std;

planet::planet(double M, double x, double y, double z, double vx, double vy, double vz)
{
    mass = M;
    position[0] = x;
    position[1] = y;
    position[2] = z;
    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;
    potential = 0.;
    kinetic = 0.;
    vec forceVector = zeros<vec>(3);
    vec previous_acceleration = zeros<vec>(3);
    vec acceleration = zeros<vec>(3);
}

planet::gForceVector(planet otherPlanet){
    // This returns the gravitational force *on* object 1 *from* object 2. Object 1
    // is pulled towards object 2.
    // Example: In sun-earth problem for the Earth orbit, the Earth is object 1 and
    // the Sun is object 2. 
    // G: gravitational constant.
    double r = norm(otherPlanet.position - this->position); // Relative distance between the objects.

    double forceStrength = (G*this->mass*otherPlanet.mass)/(r*r);   // Newton's gravitational law.

    vec forceDirection = (otherPlanet.position-this->position)/norm(otherPlanet.position-this->position);   // This vector points *from*
    // object 1 and *towards* object 2, meaning that object 1 is influenced by object 2.
    forceVector += forceStrength * forceDirection;
}