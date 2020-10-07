#include "solver_our.hpp"
#include "planet_ours.hpp"
#include <iostream>
#include <cmath>
#include "time.h"

void Solver::init(){
    vector<planet> all_planets;
    int total_planets = 0;
}

void Solver::add(planet &otherPlanet){
    all_planets.pushback(otherPlanet)
    total_planets += 1;
}

void Solver::run_velocityVerlet(double tFinal, double dt, double G){

    int N = round(tFinal/dt);   // Number of timesteps.
    int dims = 3;

    // Calculate the initial acceleration of all planets.
    for (int j=0; j<total_planets; j++){
        planet &current = all_planets[j];

        // Calculate force between current and all other planets.
        for (int k=0; k<total_planets; k++){
            // Skip if j == k.
            if (j==k){
                continue;
            }

            // This is the other planet we are comparing current to.
            planet &other = all_planets[k];
            planet::gForceVector(other);
            current.acceleration = current.forceVector / current.mass;
        }
    }
    // Loop for each time step.
    for(int i=1; i<=N-1; i++){

        // Evaluate the new position for all planets.
        for (int j=0; j<total_planets; j++){
            planet &current = all_planets[j];
            current.previous_acceleration = current.acceleration;
            current.position += dt*current.velocity + 0.5*dt*dt*current.previous_acceleration;
        }

        // Evaluate the new acceleration for all planets.
        for (int j=0; j<total_planets; j++){
            planet &current = all_planets[j];

            // Calculate the force between current and all other planets.
            for (int k=0; k<total_planets; k++){
                // Skip if j == k.
                if (j==k){
                    continue;
                }
                // This is the other planet we are comparing current to.
                planet &other = all_planets[k];
                planet::gForceVector(other);
                current.acceleration = current.forceVector / current.mass;
            }
        }

        // Evaluate the new velocity for all planets.
        for (int j=0; j<total_planets; j++){
            planet &current = all_planets[j];
            current.velocity += 0.5*dt*(current.acceleration + current.previous_acceleration);
        }
    }
}