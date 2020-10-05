#ifndef PLANET_HPP
#define PLANET_HPP

#include <cmath>
#include <armadillo>

class planet
{
public:
    // Properties
    double mass;
    double position[3];
    arma::vec position(3);
    arma::vec velocity(3);
    double potential;
    double kinetic;

    // Initializers
    planet();
    planet(double M, double x, double y, double z, double vx, double vy, double vz);

    // Functions
    double distance(planet otherPlanet);
    double GravitationalForce(planet otherPlanet, double Gconst);
    double Acceleration(planet otherPlanet, double Gconst);
    double KineticEnergy();
    double PotentialEnergy(planet &otherPlanet, double Gconst, double epsilon);

};

#endif