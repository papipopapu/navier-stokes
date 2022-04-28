#include "navierFluid.cpp"
#include "math.h"
using namespace navier;

const int scale = 1.5; const int Ni = 120 * scale + 1; const int Nj = 300*scale;
const int thick = scale * 8, midI = Ni/2, midJ = Nj - 4* thick, trueR = 10, trueWidth = 2*trueR + 1,
 lBack = 2*thick + trueWidth+1, I0 = midI+trueR+thick+1, J0 = midJ, tinyThick = thick, trueTinyLength = 5*scale, tinyLength = trueTinyLength - tinyThick + thick,
 trueLength = 60, length = trueLength + tinyThick + thick,  hJ0 = J0 - length, tinyJ0 = hJ0-thick+tinyThick;

void CI (grid &s, grid &p, grid &u, grid &v) 
{
    int tW = trueR, pW = trueR, tM = 10, pM = 0;
    line tinte(midI, midJ-thick-tW-tM, 1, false, tW);
    line presion(midI, midJ-thick-pW-pM, 1, false, pW);
    p.drawLine(presion,12000); 
    s.drawLine(tinte, 255);
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

    fluid myFluid(Ni, Nj);
    line back(I0, J0, lBack, true, thick), top(I0-lBack, hJ0, length, false, thick), bot(I0, J0-length, length, false, thick),
         mBack(I0+2*thick, J0+2*thick, lBack + 4*thick+1, true, thick),mTop(I0-2*thick-lBack, hJ0, length, false, thick), mBot(I0+2*thick, hJ0, length, false, thick)
       , tinyTop(I0-lBack+tinyLength, tinyJ0, tinyLength, true, tinyThick), tinyBot(I0, tinyJ0, tinyLength, true, tinyThick), test(midI, midJ-100, 100, false, 1);
    myFluid.drawBorders({top, back, bot, tinyTop, tinyBot, mTop, mBack, mBot});
    myFluid.boundU = 0; myFluid.boundV = 0; myFluid.boundS = 3; myFluid.boundP = 3; // free surface approx
    simulation mySimulation(myFluid);
    //mySimulation.g= 1;
    mySimulation.gradient = true;
    mySimulation.normalize = true;
    mySimulation.playSpeed = 0.5;
    mySimulation.radius = 5;
    mySimulation.draw = 0;
    mySimulation.drawMin = 0;
    mySimulation.drawMax = 0.4; // uv = 230, s = 230, p 1000, w= 0.4
    mySimulation.initialConditions = &CI;
    mySimulation.screenFactor = 6;
    //mySimulation.periodicV = &g;
    //std::cout << myFluid.solidCoordinates.is.size();
    mySimulation.start();
}