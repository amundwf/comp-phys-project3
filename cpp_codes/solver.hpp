#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "planet.hpp"
#include <fstream>
#include <armadillo>

class Solver :: public Planet{
private:
    // Parameters
    std::vector<Planet> all_planets;
    int total_planets;

public:
    //friend class Planet;
    // Functions
    void init();
    void add(Planet newPlanet);
    void run_velocityVerlet(double tFinal, double dt, double G);

};
#endif