#ifndef SOLARSYSTEM_HPP
#define SOLARSYSTEM_HPP

#include <armadillo>
#include <cmath>
#include <iostream>
#include "planet.hpp"
#include "utils.hpp"

class SolarSystem
{
public:
    SolarSystem();
    SolarSystem(arma::vec planets);

    void addPlanetToSolarSystem(planet);
    void calculateForcesAndEnergy();

    double potentialEnergy() const;
    double kineticEnergy() const;
    double totalEnergy() const;

private:
    int m_numberOfPlanets;
    arma::Col<planet> m_planets; // Vector with planets.
};

#endif