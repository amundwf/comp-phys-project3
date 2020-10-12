#include <fstream>
#include <armadillo>
#include <iostream>
#include <cmath>
#include "time.h"
#include "solver.hpp"
#include "planet.hpp"

using namespace std;
using namespace arma;

void Solver::init(){
    total_planets = 0;
}


void Solver::add(Planet newPlanet){
    // cout << this->total_planets << endl;
    //this->total_planets += 1;
    total_planets += 1;
    all_planets.push_back(newPlanet);
}
/*
void Solver::add(Planet newPlanet){
    total_planets += 1;
    all_planets.resize(total_planets);
    all_planets(total_planets-1) = newPlanet;
}
*/

int Solver::get_total_planets(){
    return total_planets;
}

std::vector<Planet> Solver::get_all_planets(){
    return all_planets;
}

void Solver::gForceVector(Planet current, Planet other, double G){
    // This returns the gravitational force *on* object 1 *from* object 2. Object 1
    // is pulled towards object 2.
    // Example: In sun-earth problem for the Earth orbit, the Earth is object 1 and
    // the Sun is object 2. 
    // G: gravitational constant.
    double r = norm(other.position - current.position); // Relative distance between the objects.

    double forceStrength = (G*current.mass*other.mass)/(r*r);   // Newton's gravitational law.

    vec forceDirection = (other.position-current.position)/norm(other.position-current.position);   // This vector points *from*
    // object 1 and *towards* object 2, meaning that object 1 is influenced by object 2.
    current.forceVector += forceStrength * forceDirection;
}

mat Solver::run_velocityVerlet(double tFinal, double dt, double G){

    int N = round(tFinal/dt);   // Number of timesteps.
    cout << "Number of timesteps: N = " << N << endl;
    //cout << "Number of planets = " << total_planets << endl;
    cout << "Number of planets = " << this->total_planets << endl;


// THE BUGS ARE SOMEWHERE WITHIN THIS COMMENT BOX

    // Set up matrix to contain all planet info.
    mat results = mat(N*(total_planets-1), 7); // Columns: t, x, y, z, vx, vy, vz
    //mat results = mat(N*total_planets, 7);
    //cout << 

    // Fill with times for all planets.
    vec tList = vec(N);
    for(int i=0; i<=N-1; i++){tList(i) = i*dt;}
    vec t_all = tList;
    for(int i=2; i < total_planets; i++){
        t_all = join_cols(t_all, tList);
    }
    results.col(0) = t_all;

    results.print("results 1:");

    // Save initial veloctiy and position to matrix.
    // Start at 1, skip the sun.

    cout << "all_planets.size(): " << all_planets.size() << endl;
    
    for (int i=0; i<=total_planets-1; i++){
        //Planet &current = all_planets[i];
        Planet current = all_planets[i];

        //Planet current = all_planets(i);
        int x = N*(i-1);
        vec currentPosition = current.getPosition();
        vec currentVelocity = current.getVelocity();

        cout << "current positon size = " << current.position.size() << endl;
        results(x, span(1,3)) = current.position.t();
        results(x, span(4,6)) = current.velocity.t();
        
    }
// END OF COMMENT BOX


    cout << "results 2\n" <<results << endl;

    // Calculate the initial acceleration of all planets.
    // start at j=1 since we don't want to update the Sun.
    for (int j=1; j<total_planets; j++){
        Planet &current = all_planets[j];

        // Calculate force between current and all other planets.
        for (int k=0; k<total_planets; k++){
            // Skip if j == k.
            if (j==k){
                continue;
            }

            // This is the other planet we are comparing current to.
            Planet &other = all_planets[k];
            Solver::gForceVector(current, other, G);
            current.acceleration = current.forceVector / current.mass;
        }
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
                Solver::gForceVector(current, other, G);
                current.acceleration = current.forceVector / current.mass;
            }
        }

        // Evaluate the new velocity for all planets.
        for (int j=0; j<total_planets; j++){
            Planet &current = all_planets[j];
            current.velocity += 0.5*dt*(current.acceleration + current.previous_acceleration);
        }
        
        // Save the new position and velocity to matrix.
        for (int j=1; j<total_planets; j++){
        Planet &current = all_planets[j];
        results(i + N*(j-1), span(1,3)) = current.position.t();
        results(i + N*(j-1), span(4,6)) = current.velocity.t();
        }
    }
    return results;
}

