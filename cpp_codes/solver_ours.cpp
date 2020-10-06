#include "solver.hpp"
#include "planet.hpp"
#include <iostream>
#include <cmath>
#include "time.h"

solver::init(){
    vector<planet> all_planets;
    int total_planets;

}

solver::add(planet &otherPlanet){
    all_planets.pushback(otherPlanet)
    total_planets += 1;
}

solver::run_velocityVerlet(double tFinal, double dt, double G){

    int N = round(tFinal/dt);   // Number of timesteps.

    // make a matrix for the acceleration to be stored.
    vec acceleration = zeros<vec>(total_planets*3);
    vec previous_acceleration = zeros<vec>(total_planets*3);

    // Loop for each time step.
    for(int i=1; i<=N-1; i++){
        // Loop for each planet.
        for (int j=0; j<total_planets; j++){
            // This is the current planet whose position and velocity we 
            // will be updating. 
            planet &current = all_planets[j];

            // Calculate force between current and all other planets.
            // not sure about the j+1, it means we don't find the force of
            // it on itself. but then it misses the planets before it.
            for (int k = j+1; k<total_planets; k++){
                // This is the other planet we are comparing current to.
                planet &other = all_planets[k];
                planet::gForceVector(other);

            previous_acceleration = acceleration();
            previous_acceleration = current.ForceVector/current.mass; 
            
            previous_acceleration = acceleration;

            // Evaluate the current position, acceleration and velocity.
            position += dt*velocity + 0.5*dt*dt*previous_acceleration;
            acceleration = gForceVector(G, m_E, m_S, position, sunPosition)/m_E;
            velocity += 0.5*dt*(acceleration + previous_acceleration);

            }
        }
    }
}