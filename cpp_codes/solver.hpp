#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "planet.hpp"
#include <fstream>
#include <armadillo>

class Solver{
private:
    // Allow Solver access to member function of planet.
    friend class Planet;

    // Parameters
    std::vector<Planet> all_planets;
    int total_planets;

public:
    // Functions
    void Solver::init();
    void Solver::add(Planet newPlanet);
    void Solver::run_velocityVerlet(double tFinal, double dt, double G);

};
#endif