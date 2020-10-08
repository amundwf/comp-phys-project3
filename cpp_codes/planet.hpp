#ifndef PLANET_HPP
#define PLANET_HPP

#include <cmath>
#include <armadillo>

class Planet{
public:
    // Parameters
    arma::vec position;
    arma::vec velocity;
    arma::vec forceVector;
    arma::vec previous_acceleration;
    arma::vec acceleration;
    double mass;
    double potential;
    double kinetic;

    // Functions
    void planet(double M, arma::vec initialPosition, arma::vec initialVelocity);
    void gForceVec(Planet otherPlanet);
};
#endif