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
    arma::mat angMomentum_energy_mat;
    int revolution;
    
    
public:
    arma::mat perihelion_mat_solver;
    // Functions
    void init(int N);
    void add(Planet newPlanet);
    arma::mat run_velocityVerletBeta(double tFinal, double dt, double beta, double G);
    arma::mat run_velocityVerletForceType(int functionNum, double tFinal, double dt, double G);
    arma::mat get_angMomentum_energy_mat();
    int get_total_planets();
    std::vector<Planet> get_all_planets();
    void totalEnergySystem(int i, double G);
    void totalAngularMomentumSystem(int i);
    arma::vec eval_perihelion(Planet &current);
};
#endif