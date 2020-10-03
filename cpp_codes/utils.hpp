#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <cmath>
#include <armadillo>


// This file contains various functions, for running different parts of the project.
arma::vec gForceVector(double G, double mass1, double mass2, arma::vec pos1, arma::vec pos2);

void run_forwardEuler();

void run_velocityVerlet();

#endif 
