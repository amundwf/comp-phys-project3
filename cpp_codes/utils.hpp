#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <cmath>
#include <armadillo>


// This file contains various functions, for running different parts of the project.
arma::vec gForceVector(double G, double mass1, double mass2, arma::vec pos1, arma::vec pos2);

void writeResultsToFile(arma::mat results, std::string filename, std::string directory);

arma::mat run_forwardEuler(double tFinal, double dt, double G);

arma::mat run_velocityVerlet(double tFinal, double dt, double G);

void task_3a_forwardEuler(double G);

//void task_3a_velocityVerlet(double G);

#endif 
