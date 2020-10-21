#include <fstream>
#include <armadillo>
#include <iostream>
#include <cmath>
#include "time.h"
#include "solver.hpp"
#include "planet.hpp"
#include "utils.hpp"

using namespace std;
using namespace arma;

void Solver::init(int N){
    total_planets = 0;
    angMomentum_energy_mat = mat(N,4);
    // angMomentum_energy_mat: 4 columns: total kinetic energy, 
    // total potential energy, total mechanical energy, total angular
    // momentum.
}

void Solver::add(Planet newPlanet){
    // cout << this->total_planets << endl;
    //this->total_planets += 1;
    total_planets += 1;
    all_planets.push_back(newPlanet);
}

int Solver::get_total_planets(){
    return total_planets;
}

std::vector<Planet> Solver::get_all_planets(){
    return all_planets;
}

void Solver::totalAngularMomentumSystem(int i){
    // Calculate the total angular momentum of the system.

    double totalAngMomentum = 0.0;
    for (int j=1; j <= total_planets-1; j++){
        Planet current = all_planets[j];
        totalAngMomentum += current.angularMomentum();
    }
    angMomentum_energy_mat(i, 3) = totalAngMomentum;
}

void Solver::totalEnergySystem(int i, double G){
    // Calculate the total energy of the system.

    // Initialise totalK to 0.
    double totalKinetic = 0.0;
    double totalPotential = 0.0;
    for (int j=1; j <= total_planets-1; j++){
        Planet current = all_planets[j];

        totalKinetic += current.kineticEnergy();

        for (int k=0; k <= total_planets-1; k++){
            // Skip if j == k.
            if (j==k){
                continue;
            }
            // This is the other planet we are comparing current to.
            Planet other = all_planets[k];
            totalPotential += potentialEnergy(current, other, G);
        }
    }
    double totalEnergy = (totalKinetic + totalPotential);
    angMomentum_energy_mat(i, 0) = totalKinetic;
    angMomentum_energy_mat(i, 1) = totalPotential;
    angMomentum_energy_mat(i, 2) = totalEnergy;
}

mat Solver::get_angMomentum_energy_mat(){
    return angMomentum_energy_mat;
}

mat Solver::run_velocityVerlet(double tFinal, double dt, double G){
    int N = round(tFinal/dt);   // Number of timesteps.
    //cout << "Number of timesteps: N = " << N << endl;
    //cout << "Number of planets = " << this->total_planets << endl;

    // Set up matrix to contain all planet info.
    mat results = mat(N*(total_planets-1), 7);

    // Fill matrix with times for all planets.
    vec tList = vec(N);
    for(int i=0; i<=N-1; i++){tList(i) = i*dt;}
    vec t_all = tList;
    for(int i=2; i < total_planets; i++){
        t_all = join_cols(t_all, tList);
    }
    results.col(0) = t_all;

    // Calculate the initial acceleration of all planets.
    // start at j=1 since we don't want to update the Sun.
    for (int j=1; j <= total_planets-1; j++){
        Planet &current = all_planets[j];

        // Calculate force between current and all other planets.
        vec totalForce = vec("0 0 0");  // Instantiate it as zero force to start
        // with and then add the forces from the sun and all planets.
        vec force;

        for (int k=0; k <= total_planets-1; k++){
            // Skip if j == k.
            if (j==k){
                continue;
            }

            // This is the other planet we are comparing current to.
            Planet &other = all_planets[k];
            // Get the force from the other planet on the current planet.
            force = gForceVectorPlanet(current, other, G);
            totalForce += force;
        }
        current.forceVector = totalForce;
        current.acceleration = current.forceVector / current.mass;

        /*
        // Print some planet info.
        cout << "\nPlanet number: " << j << endl;
        cout <<"Intial acceleration: " << current.acceleration.t() << endl;
        cout <<"Intial forceVector: " << current.forceVector.t() << endl;
        cout <<"current.mass: " << current.mass << endl;
        */
    }

    // Calculate initial total energy.
    totalEnergySystem(0, G);
    totalAngularMomentumSystem(0);

    // Save initial velocity and position and total energy to matrix.
    // Start at 1, skip the sun.
    for (int i=1; i<total_planets; i++){
        //Planet &current = all_planets[i];
        Planet current = all_planets[i];

        //vec currentPosition = current.getPosition();
        //vec currentVelocity = current.getVelocity();
        int x = N*(i-1);
        results(x, span(1,3)) = current.position.t();
        results(x, span(4,6)) = current.velocity.t();
    }

    // Loop for each time step.
    for(int i=1; i<=N-1; i++){

        // Evaluate the new position for all planets.
        // start at j=1 since we don't want to update the Sun.
        for (int j=1; j<total_planets; j++){
            Planet &current = all_planets[j];
            current.previous_acceleration = current.acceleration;
            current.position += dt*current.velocity + 0.5*dt*dt*current.previous_acceleration;
        }

        // Evaluate the new acceleration for all planets.
        vec totalForce = vec("0 0 0"); // Initialize sum of forces as zero.
        vec force;
        for (int j=1; j<total_planets; j++){
            Planet &current = all_planets[j];

            // Calculate the force between current and all other planets.
            for (int k=0; k<total_planets; k++){
                // Skip if j == k.
                if (j==k){
                    continue;
                }
                // This is the other planet we are comparing current to.
                Planet &other = all_planets[k];
                force = gForceVectorPlanet(current, other, G);
                totalForce += force;
            }
            current.forceVector = totalForce;
            current.acceleration = current.forceVector / current.mass;
        }

        // Evaluate the new velocity for all planets.
        for (int j=1; j<total_planets; j++){
            Planet &current = all_planets[j];
            current.velocity += 0.5*dt*(current.acceleration + current.previous_acceleration);
        }
        
        // Save the new position and velocity to matrix.
        for (int j=1; j<total_planets; j++){
            Planet &current = all_planets[j];
            results(i + N*(j-1), span(1,3)) = current.position.t();
            results(i + N*(j-1), span(4,6)) = current.velocity.t();
        }

        // Print the total energy to results.
        totalEnergySystem(i, G);
        totalAngularMomentumSystem(i);
    }
    return results;
}

mat Solver::run_velocityVerletBeta(double tFinal, double dt, double beta, double G){
    int N = round(tFinal/dt);   // Number of timesteps.
    //cout << "Number of timesteps: N = " << N << endl;
    //cout << "Number of planets = " << total_planets << endl;

    // Set up matrix to contain all planet info.
    mat results = mat(N*(total_planets-1), 7);

    // Fill matrix with times for all planets.
    vec tList = vec(N);
    for(int i=0; i<=N-1; i++){tList(i) = i*dt;}
    vec t_all = tList;
    for(int i=2; i < total_planets; i++){
        t_all = join_cols(t_all, tList);
    }
    results.col(0) = t_all;

    // Calculate the initial acceleration of all planets.
    // start at j=1 since we don't want to update the Sun.
    for (int j=1; j <= total_planets-1; j++){
        Planet &current = all_planets[j];

        // Calculate force between current and all other planets.
        vec totalForce = vec("0 0 0");  // Instantiate the total force as zero
        // to start with and then add the forces from the sun and all planets.
        vec force;
        // This loop includes the Sun.
        for (int k=0; k <= total_planets-1; k++){
            // Skip if j == k.
            if (j==k){
                continue;
            }

            // This is the other planet we are comparing current to.
            Planet &other = all_planets[k];
            // Get the force from the other planet on the current planet.
            force = gForcePlanetBeta(current, other, beta, G);
            totalForce += force;
        }
        current.forceVector = totalForce;
        current.acceleration = current.forceVector / current.mass;

        /*
        // Print some planet info.
        cout << "\nPlanet number: " << j << endl;
        cout <<"Intial acceleration: " << current.acceleration.t() << endl;
        cout <<"Intial forceVector: " << current.forceVector.t() << endl;
        cout <<"current.mass: " << current.mass << endl;
        */
    }
    // Calculate initial total energy.
    totalEnergySystem(0, G);
    totalAngularMomentumSystem(0);

    // Save initial velocity and position and total energy to matrix.
    // Start at 1, skip the sun.
    for (int i=1; i<total_planets; i++){
        //Planet &current = all_planets[i];
        Planet current = all_planets[i];

        //vec currentPosition = current.getPosition();
        //vec currentVelocity = current.getVelocity();
        int x = N*(i-1);
        results(x, span(1,3)) = current.position.t();
        results(x, span(4,6)) = current.velocity.t();
    }

    // Loop through all timesteps:
    for(int i=1; i<=N-1; i++){
        // Evaluate the new position for all planets.
        // start at j=1 since we don't want to update the Sun.
        for (int j=1; j<total_planets; j++){
            Planet &current = all_planets[j];
            current.previous_acceleration = current.acceleration;
            current.position += dt*current.velocity + 0.5*dt*dt*current.previous_acceleration;
        }

        // Evaluate the new acceleration for all planets.
        vec totalForce = vec("0 0 0"); // Initialize sum of forces as zero.
        vec force;
        for (int j=1; j<total_planets; j++){
            Planet &current = all_planets[j];

            // Calculate the force between current and all other planets.
            for (int k=0; k<total_planets; k++){
                // Skip if j == k.
                if (j==k){
                    continue;
                }
                // This is the other planet we are comparing current to.
                Planet &other = all_planets[k];
                force = gForcePlanetBeta(current, other, beta, G);
                totalForce += force;
            }
            current.forceVector = totalForce;
            current.acceleration = current.forceVector / current.mass;
        }

        // Evaluate the new velocity for all planets.
        for (int j=1; j<total_planets; j++){
            Planet &current = all_planets[j];
            current.velocity += 0.5*dt*(current.acceleration + current.previous_acceleration);
        }
        
        // Save the new position and velocity to matrix.
        for (int j=1; j<total_planets; j++){
            Planet &current = all_planets[j];
            results(i + N*(j-1), span(1,3)) = current.position.t();
            results(i + N*(j-1), span(4,6)) = current.velocity.t();
        }

        // Print the total energy to results.
        totalEnergySystem(i, G);
        totalAngularMomentumSystem(i);
    }
    return results;
}
