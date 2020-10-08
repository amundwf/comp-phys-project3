#include "solarsystem.hpp"

using namespace std;
using namespace arma;

SolarSystem::SolarSystem() :
    m_kineticEnergy(0),
    m_potentialEnergy(0)
{}

/*
void calculateForcesAndEnergy(){
    // Calculates the forces and energy for the entire
    // solar system. (Why do we need this?)
    m_kineticEnergy = 0;
    m_potentialEnergy = 0;
}
*/

/*
void addPlanetToSolarSystem(planet p){
    m_planets.append(planet)
}
*/