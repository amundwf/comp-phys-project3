#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <fstream>
#include <armadillo>
#include "planet.hpp"

class Solver{
private:

    int total_planets = 0;
    std::vector<Planet> all_planets;
    
public:

    // Functions
    void init();
    void add(Planet newPlanet);
    void gForceVector(Planet current, Planet other, double G);
    arma::mat run_velocityVerlet(double tFinal, double dt, double G);

};
#endif