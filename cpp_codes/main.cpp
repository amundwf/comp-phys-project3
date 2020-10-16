#include <iostream>
#include <cmath>
#include "utils.hpp"

using namespace std;
using namespace arma;

int main(){
    double G = 4*M_PI*M_PI; // The gravitational constant when having replaced the units
    // kg, meters and seconds by the units solar masses, astronomical units and years, respectively.
    
    //task_3a_forwardEuler(G);

    //task_3a_velocityVerlet(G);

    task_3b_velocityVerlet(G);
    
    //double initialSpeed_kmPerSec = 40;
    //task_3f_escape_velocity(initialSpeed_kmPerSec, G);

    //task_3g_three_body(G);

    return 0;
}