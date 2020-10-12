#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <fstream>
#include <armadillo>
#include <iostream>
#include <cmath>
#include "time.h"
#include "solver.hpp"
#include "planet.hpp"


class Solver{
private:
    int total_planets;
    std::vector<Planet> all_planets;
    //arma::Col<Planet> all_planets;
    
public:

    // Functions
    void init();
    void add(Planet newPlanet);
    void gForceVector(Planet current, Planet other, double G);
    arma::mat run_velocityVerlet(double tFinal, double dt, double G);
    int get_total_planets();
    std::vector<Planet> get_all_planets();

};
#endif