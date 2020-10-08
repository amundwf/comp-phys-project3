#ifndef PLANET_HPP
#define PLANET_HPP

#include <cmath>
#include <armadillo>

class planet{
public:
    // Parameters
    arma::vec position(3);
    arma::vec velocity(3);
    arma::vec forceVector(3);
    arma::vec previous_acceleration(3);
    arma::vec acceleration(3);
    double mass;
    double potential;
    double kinetic;

    // Functions
    planet(double M, double x, double y, double z, double vx, double vy, double vz);
    gForceVec(planet otherPlanet);
};
#endif