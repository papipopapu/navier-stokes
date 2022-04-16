#include <cstring>
#include <iostream>
#include <algorithm>
#include <SFML/Graphics.hpp>


namespace navier
{
class line
{
    public:
    int i0, j0, length, thickness;
    bool isVertical;
    line ();
    line (int i0_, int j0_, int length_, bool isVertical_, int thickness);
};


class grid
{
    private:
    size_t Ni, Nj;
    std::vector<float> V; 
    public:
    grid ();
    grid (size_t Ni_, size_t Nj_);
    void drawLine(line Line, float val); 
    void fill (size_t Ni_, size_t Nj_);
    std::vector<float> & getV();
    const std::vector<float>  & getV()  const;
    float & operator()(size_t i, size_t j);
    const float operator()(size_t i, size_t j) const;

    const size_t & getNi() const;
    const size_t & getNj() const;
    size_t & getNi();
    size_t & getNj();
};

class coordinates
{
    public:
    std::vector<int> is, js;
    coordinates();
    coordinates(grid &map);
    void fill(grid &map);
};
void pass(grid &s, double t, double dt) {}
void pass(grid &s, grid &p, grid &u, grid &v) {}
void paintPixels(grid &s, grid &solidMap, std::vector<sf::Uint8> &pixels, float r, float g, float b, int a, float max, float min); // to's end's not included
void addSource(int radius, grid &S, int i, int j, int b, float val, float maxVal, coordinates &solidCoordinates, grid &solidMap);
void advect(int b, int solidB, float dt, grid &s, grid &s0, grid &u, grid &v, coordinates &solidCoordinates, grid &solidMap);
void diffuse(int b, int solidB, float dt, float diff, grid &s, grid &s0, coordinates &solidCoordinates, grid &solidMap);
void project(int bu, int bv, int bp, grid &u, grid &v, grid &p, grid &div, coordinates &solidCoordinates, grid &solidMap);
void dissipate(grid &s, double dt, float disip);
void boundaries(int b, grid &x);
void lineBoundaries(int b, grid &x, coordinates &solidCoordinates, grid &solidMap);
void swapV(grid& u1, grid& u2);

class fluid
{
    private:
    int Ni, Nj;
    float t = 0; // normalized dh = 1
    grid u, v, u0, v0, s, s0, p, solidMap;
    coordinates solidCoordinates;

    public:
    int boundU = 1, boundV = 2, boundP = 0, boundS = 0;
    float dt, visc = 1, diff = 1, disip = 0.005; 
    fluid(int Ni_, int Nj_);// Ni, Nj refer to the dimensions including solid cells
    void drawBorders(std::initializer_list<line> Lines);
    void addS(int i, int j, int radius, float val, float maxVal);
    void addP(int i, int j, int radius, float val, float maxVal);
    void addV(int i, int j, int radius, float uVal, float vVal, float maxValU, float maxValV);
    void vStep();
    void sStep();
    void step(float dt_);

    friend class simulation;

};
class simulation
{
    private:
    fluid Fluid;

    size_t Ni, Nj;
    public:
    void (*periodicU) (grid &s, double t, double dt) = &pass,
         (*periodicV) (grid &s, double t, double dt) = &pass,
         (*periodicS) (grid &s, double t, double dt) = &pass,
         (*periodicP) (grid &s, double t, double dt) = &pass;

    void (*initialConditions) (grid &s, grid &p, grid &u, grid &v) = &pass;
    float mouseS = 2000, mouseV = 10, mouseP = 100000; // increment of addition
    float maxAddS = 2000, maxAddU = 100, maxAddV = 100, maxAddP = 100000;
    float r = 0, g = 0, b = 1;
    float screenFactor = 1;
    float playSpeed = 0.2;
    bool normalize = true;
    uint8_t a = 255; // colors
    uint8_t draw = 0, radius = 3;// 0 for s, 1 for u, 2 for v, 3 for p
    uint8_t fps = 33;

    simulation (fluid Fluid);
    void start();
};
}