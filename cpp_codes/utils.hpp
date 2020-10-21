#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <cmath>
#include <armadillo>
#include "solver.hpp" 
#include "planet.hpp"

// This file contains various functions, for running different parts of the project.
arma::vec gForceVector(double G, double mass1, double mass2, arma::vec pos1, arma::vec pos2);
arma::vec gForceVectorPlanet(Planet planet1, Planet planet2, double G);
arma::vec gForcePlanetBeta(Planet planet1, Planet planet2, double beta, double G);

double potentialEnergy(Planet current, Planet other, double G);

void writeMatrixToFile(arma::mat results, std::string filename, std::string directory);
void writeGeneralMatrixToCSV(arma::mat results, arma::field<std::string> columnLabels, std::string filename, std::string directory);

arma::vec initial_earth_velocity(arma::vec initialPosition);

double get_earth_mass();

// Convertors between velocity units:
double kmPerSec_to_auPerYear(double speed_kmPerSec);
double auPerYear_to_kmPerSec(double speed_auPerYear);

arma::mat run_forwardEuler(double tFinal, double dt, double G);
arma::mat forwardEuler(double tFinal, double dt, double m_SI, arma::vec initialPosition, arma::vec initialVelocity, double G);

arma::mat run_velocityVerlet(double tFinal, double dt, double G);
arma::mat velocityVerlet(double tFinal, double dt, double m_SI, arma::vec initialPosition, arma::vec initialVelocity, double G);

void task_3a_forwardEuler(double G);

void task_3a_velocityVerlet(double G);

void task_3b_velocityVerlet(double G);

void task_3e_force(double G);

void task_3f_escape_velocity(double initialSpeed_kmPerSec, double G);

void task_3g_three_body(double G);

//void run_solarSystem();

#endif 
