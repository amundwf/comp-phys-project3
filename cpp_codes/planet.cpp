#include "planet.hpp"
#include <armadillo>

using namespace arma;
using namespace std;

void Planet::init(double M, vec initialPosition, vec initialVelocity){
    mass = M;
    potential = 0.;
    kinetic = 0.;
    position = initialPosition;
    velocity = initialVelocity;
    forceVector = zeros<vec>(3);
    previous_acceleration = zeros<vec>(3);
    acceleration = zeros<vec>(3);
}

vec Planet::getPosition(){
    return position;
}

vec Planet::getVelocity(){
    return velocity;
}

double Planet::kineticEnergy(){
    return 0.5*mass*dot(velocity.t(),velocity);
}

double Planet::angularMomentum(){
    // This function returns the absolute value of the angular momentum
    // of the planet.
    return mass*norm(cross(position, velocity));
}

