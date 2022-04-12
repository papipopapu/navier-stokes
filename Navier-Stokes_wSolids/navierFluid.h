#include <cstring>
#include <iostream>
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
    size_t getNi() const;
    size_t getNj() const;
};

class coordinates
{
    public:
    std::vector<int> is, js;
    coordinates();
    coordinates(grid &map);
    void fill(grid &map);
};
void paintPixels(grid &s, grid &solidMap, std::vector<sf::Uint8> &pixels, float r, float g, float b, int a); // to's end's not included
void addSource(int radius, grid &S, int i, int j, int b, float val, coordinates &solidCoordinates, grid &solidMap);
void advect(int b, float dt, grid &s, grid &s0, grid &u, grid &v, coordinates &solidCoordinates, grid &solidMap);
void diffuse(int b, float dt, float diff, grid &s, grid &s0, coordinates &solidCoordinates, grid &solidMap);
void project(grid &u, grid &v, grid &p, grid &div, coordinates &solidCoordinates, grid &solidMap);
void boundaries(int b, grid &x);
void lineBoundaries(int b, grid &x, coordinates &solidCoordinates, grid &solidMap);
void swapV(grid& u1, grid& u2);

class fluid
{
    public:
    int Ni, Nj;
    float t = 0, dt, visc = 1, diff = 1; // normalized dh = 1
    grid u, v, u0, v0, s, s0, p, solidMap;
    coordinates solidCoordinates;

    fluid(int Ni_, int Nj_);// Ni, Nj refer to the dimensions including solid cells
    void drawBorders(std::initializer_list<line> Lines);
    void addS(int i, int j, int radius, float val);
    void addP(int i, int j, int radius, float val);
    void addV(int i, int j, int radius, float uVal, float vVal);
    void vStep();
    void sStep();
    void step(float dt_);
    size_t getNi() const;
    size_t getNj() const;

};
class simulation
{
    private:
    fluid Fluid;

    size_t WIDTH, HEIGHT, Ni, Nj;
    sf::Texture texture;
    sf::Sprite sprite;
    sf::Image image;
    sf::RenderWindow window;

    public:
    float mouseS = 500, mouseV = 10, mouseP = 50000; // increment of addition
    float r = 0, g = 0, b = 1;
    uint8_t a = 255; // colors
    uint8_t draw = 0, radius = 3;// 0 for s, 1 for u, 2 for v, 3 for p
    uint8_t fps = 33;

    simulation (fluid Fluid);
    void start();
};
}