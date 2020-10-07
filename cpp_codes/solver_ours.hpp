#ifndef SOLVER_OURS_HPP
#define SOLVER_OURS_HPP

#include "planet_ours.hpp"
#include <fstream>
#include <armadillo>

class Solver{
private:
    // Allow Solver access to member function of planet.
    friend class planet;
    // Parameters
    vector<planet> all_planets;
    int total_planets;

public:
    // Functions
    void Solver::init();
    void Solver::add(planet otherPlanet);
    void Solver::run_veloctityVerlet(double tFinal, double dt, double G)

};
#endif