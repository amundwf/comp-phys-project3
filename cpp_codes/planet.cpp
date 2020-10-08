#include "planet.hpp"
#include <armadillo>

using namespace arma;
using namespace std;

Planet::init(double M, vec initialPosition, vec initialVelocity){
    double mass = M;
    double potential = 0.;
    double kinetic = 0.;
    vec position = initialPosition;
    vec velocity = initialVelocity;
    vec forceVector = zeros<vec>(3);
    vec previous_acceleration = zeros<vec>(3);
    vec acceleration = zeros<vec>(3);
}

Planet::gForceVector(Planet otherPlanet){
    // This returns the gravitational force *on* object 1 *from* object 2. Object 1
    // is pulled towards object 2.
    // Example: In sun-earth problem for the Earth orbit, the Earth is object 1 and
    // the Sun is object 2. 
    // G: gravitational constant.
    double r = norm(otherPlanet.position - position); // Relative distance between the objects.

    double forceStrength = (G*mass*otherPlanet.mass)/(r*r);   // Newton's gravitational law.

    vec forceDirection = (otherPlanet.position-position)/norm(otherPlanet.position-position);   // This vector points *from*
    // object 1 and *towards* object 2, meaning that object 1 is influenced by object 2.
    forceVector += forceStrength * forceDirection;
}