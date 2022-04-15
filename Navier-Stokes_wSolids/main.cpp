#include "navierFluid.cpp"
#include "math.h"
using namespace navier;


void CI (grid &s, grid &p, grid &u, grid &v) 
{
    int length = 70, width = 40, i0 = 80, j0 = 280, dlength = 10, thick = 8, extraTinte = 5, mueveTinte = 15;
    int zone = 10;
    line tinte(i0-(int)width/2., j0-zone-mueveTinte, 1, false, zone+extraTinte);
    line presion(i0-(int)width/2., j0-zone, 1, false, zone);
    p.drawLine(presion, 15000); s.drawLine(tinte, 255);
}
void g(grid &v, double t, double dt)
{
    for (int i=0;  i<v.getNi()*v.getNj(); i++)
    {
        v.getV()[i] -= 9.8 * dt;
    }
}
int main()
{
    int length = 70, width = 40, i0 = 80, j0 = 280, dlength = 10, thick = 8, tiny = 5;
    fluid myFluid(120, 300);
    line back(i0, j0, width+1, true, thick), top(i0-width, j0-length, length, false, thick), bot(i0, j0-length, length, false, thick)
       , tinyTop(i0+tiny-width, j0-length, tiny, true, thick), tinyBot(i0, j0-length, tiny, true, thick); 
    myFluid.drawBorders({top, back, bot, tinyTop, tinyBot});
    simulation mySimulation(myFluid);
    mySimulation.g = 1;
    mySimulation.radius = 5;
    mySimulation.draw = 3;
    mySimulation.initialConditions = &CI;
    mySimulation.screenFactor = 5;
    mySimulation.periodicV = &g;
    //std::cout << myFluid.solidCoordinates.is.size();
    mySimulation.start();
}