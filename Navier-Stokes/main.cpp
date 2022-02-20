#include "navierFluid.cpp"
using namespace navier;
int main()
{
    fluid myFluid(202, 202);
    simulation mySimulation(myFluid);
    mySimulation.g = 1;
    mySimulation.b = 0;
    mySimulation.radius = 5;
    mySimulation.draw = 0;
    mySimulation.start();
}