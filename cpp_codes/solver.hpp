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
    arma::mat energyMatrix;
public:
    // Functions
    void init(int N);
    void add(Planet newPlanet);
    arma::mat run_velocityVerlet(double tFinal, double dt, double G);
    arma::mat get_energy_matrix();
    int get_total_planets();
    std::vector<Planet> get_all_planets();
    void totalEnergySystem(int i, double G);
};
#endif