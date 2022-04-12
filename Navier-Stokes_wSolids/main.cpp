#include "navierFluid.cpp"
using namespace navier;
int main()
{
    int length = 70, width = 40, i0 = 120, j0 = 180, dlength = 10, thick = 8;
    fluid myFluid(202, 202);
    line back(i0, j0, width, true, thick), top(i0-width, j0-length, length, false, thick), bot(i0, j0-length, length, false, thick);
    myFluid.drawBorders({top, back, bot});
    simulation mySimulation(myFluid);
    mySimulation.g = 1;
    mySimulation.b = 0;
    mySimulation.radius = 5;
    mySimulation.draw = 0;
    //std::cout << myFluid.solidCoordinates.is.size();
    mySimulation.start();
}