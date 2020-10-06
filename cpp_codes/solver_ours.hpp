#ifndef SOLVER_OURS_HPP
#define SOLVER_OURS_HPP

#include "planet_ours.hpp"
#include <fstream>
#include <armadillo>

class solver
{
public:
    friend class planet;

    vector<planet> all_planets;
    int total_planets;
};

#endif