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

    for(int i=1; i<=N-1; i++){
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

            
            }
        }
    }
}