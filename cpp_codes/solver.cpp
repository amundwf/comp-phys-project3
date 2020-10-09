#include "solver.hpp"

#include <iostream>
#include <cmath>
#include "time.h"

using namespace std;
using namespace arma;

void Solver::init(){
    vector<Planet> all_planets;
    int total_planets = 0;
}

void Solver::add(Planet &newPlanet){
    all_planets.pushback(newPlanet)
    total_planets += 1;
}

void Solver::run_velocityVerlet(double tFinal, double dt, double G){

    int N = round(tFinal/dt);   // Number of timesteps.
    int dims = 3;

    // Calculate the initial acceleration of all planets.
    for (int j=0; j<total_planets; j++){
        Planet &current = all_planets[j];

        // Calculate force between current and all other planets.
        for (int k=0; k<total_planets; k++){
            // Skip if j == k.
            if (j==k){
                continue;
            }

            // This is the other planet we are comparing current to.
            Planet other = all_planets[k];
            current.gForceVector(other);
            current.acceleration = current.forceVector / current.mass;
        }
    }
    // Loop for each time step.
    for(int i=1; i<=N-1; i++){

        // Evaluate the new position for all planets.
        for (int j=0; j<total_planets; j++){
            Planet &current = all_planets[j];
            current.previous_acceleration = current.acceleration;
            current.position += dt*current.velocity + 0.5*dt*dt*current.previous_acceleration;
        }

        // Evaluate the new acceleration for all planets.
        for (int j=0; j<total_planets; j++){
            Planet &current = all_planets[j];

            // Calculate the force between current and all other planets.
            for (int k=0; k<total_planets; k++){
                // Skip if j == k.
                if (j==k){
                    continue;
                }
                // This is the other planet we are comparing current to.
                Planet other = all_planets[k];
                current.gForceVector(other);
                current.acceleration = current.forceVector / current.mass;
            }
        }

        // Evaluate the new velocity for all planets.
        for (int j=0; j<total_planets; j++){
            Planet &current = all_planets[j];
            current.velocity += 0.5*dt*(current.acceleration + current.previous_acceleration);
        }
        
        cout << all_planets(0).position << all_planets(0).velocity << endl;
        cout << all_planets(1).position << all_planets(1).velocity << endl;
        
    }
}