#include <iostream>
#include <cmath>
//#include "planet.hpp"
//#include "solver.hpp"
#include "utils.hpp"

using namespace std;
using namespace arma;

int main(){
    double G = 6.67e-11;

    //run_forwardEuler();
    vec force1 = gForceVector(G, 1/double(1000000), 1, vec("1 0 0"), vec("0 0 0"));
    force1.print("force1:");

    //run_velocityVerlet();

    return 0;
}