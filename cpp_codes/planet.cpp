#include "planet.hpp"
#include <armadillo>

using namespace arma;
using namespace std;

void Planet::init(double M, vec initialPosition, vec initialVelocity){
    double mass = M;
    double potential = 0.;
    double kinetic = 0.;
    vec position = initialPosition;
    cout << "intial postion size = " << position.size() <<endl;
    vec velocity = initialVelocity;
    vec forceVector = zeros<vec>(3);
    vec previous_acceleration = zeros<vec>(3);
    vec acceleration = zeros<vec>(3);
}

