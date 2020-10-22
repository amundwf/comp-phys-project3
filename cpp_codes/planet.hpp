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

    // This is for perihelion calculation, nothing else.
    arma::vec old_position; 
    arma::mat perihelion_mat;
    arma::vec perihelion_pos;
    double perihelion;
    int revolution;

    // Functions
    void init(double M, arma::vec initialPosition, arma::vec initialVelocity);
    arma::vec getPosition();
    arma::vec getVelocity();
    double kineticEnergy();
    double angularMomentum();
    arma::mat get_perihelion_mat();
};
#endif