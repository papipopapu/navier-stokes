#include "navierFluid.h"

navier::line::line() {}
navier::line::line(int i0_, int j0_, int length_, bool isVertical_, int thickness_) 
: i0(i0_), j0(j0_), length(length_), isVertical(isVertical_), thickness(thickness_) {}
navier::coordinates::coordinates() {}
navier::coordinates::coordinates(grid &map)
{
    int i,j;
    for (i=0; i++; i<map.getNi())
    {
        for (j=0; j++; j<map.getNj())
        {
            if (map(i, j) != 0)
            {
                is.push_back(i); js.push_back(j);
            }
        }
    }
}
void navier::coordinates::fill(grid &map)
{
    int i,j;
    for (i=0; i<map.getNi(); i++)
    {
        for (j=0; j<map.getNj(); j++)
        {
            if (map(i, j) != 0)
            {
                is.push_back(i); js.push_back(j);
                
            }
            
        }
    }
}
navier::grid::grid () {}
navier::grid::grid (size_t Ni_, size_t Nj_) : Nj(Nj_), Ni(Ni_)
{
    V.resize(Nj*Ni);
}
void navier::grid::drawLine (line Line, float val) // line vector in direction (ux, uy)
{
    int i = Line.i0, j = Line.j0, l = Line.thickness;
    for(int d=0; d<Line.length; d++)
    {
        (*this)(i, j) = val;

        for (int ip=-l; ip<=l; ip++)
        {
            for (int jp=-l; jp<=l; jp++)
            {
                (*this)(i+ip, j+jp) = val;
            }
        }

        if (Line.isVertical) {i--;}
        else {j++;}
    }
}
void navier::grid::fill (size_t Ni_, size_t Nj_) 
{
    Nj = Nj_; Ni = Ni_;
    V.resize(Nj*Ni);
}
std::vector<float> & navier::grid::getV() 
{
    return V;
}
const std::vector<float> & navier::grid::getV()  const
{
    return V;
}
float & navier::grid::operator()(size_t i, size_t j) 
{
    return V[i*Nj + j];
}
const float navier::grid::operator()(size_t i, size_t j) const
{
    return V[i*Nj + j];
}
const size_t & navier::grid::getNi() const {return Ni;}
const size_t & navier::grid::getNj() const {return Nj;}
size_t & navier::grid::getNi()  {return Ni;}
size_t & navier::grid::getNj()  {return Nj;}


void navier::swapV(grid& u1, grid& u2) { std::swap(u1.getV() , u2.getV());}
navier::fluid::fluid(int Ni_, int Nj_) : Ni(Ni_), Nj(Nj_) // Ni, Nj refer to the dimensions including solid cells
{       
    u.fill(Ni, Nj); u0.fill(Ni, Nj); v.fill(Ni, Nj); v0.fill(Ni, Nj);
    s.fill(Ni, Nj); s0.fill(Ni, Nj); p.fill(Ni, Nj); solidMap.fill(Ni, Nj);
}
void navier::fluid::drawBorders(std::initializer_list<line> Lines)
{
    for (line L : Lines)
    {
        solidMap.drawLine(L, 1);
    }
    solidCoordinates.fill(solidMap);
}
void navier::fluid::addS(int i, int j, int radius, float val)
{
    addSource(radius, s, i, j, 0, val, solidCoordinates, solidMap);
}
void navier::fluid::addP(int i, int j, int radius, float val)
{
    addSource(radius, p, i, j, 0, val, solidCoordinates, solidMap);
}
void navier::fluid::addV(int i, int j, int radius, float uVal, float vVal)
{
    addSource(radius, u, i, j, 1, uVal, solidCoordinates, solidMap);
    addSource(radius, v, i, j, 2, vVal, solidCoordinates, solidMap);
}
void navier::fluid::vStep()
{
    swapV(u, u0); swapV(v, v0);
    diffuse(2, dt, visc, v, v0, solidCoordinates, solidMap); diffuse(1, dt, visc, u, u0, solidCoordinates, solidMap);
    project(u, v, p, v0, solidCoordinates, solidMap); // MOTHERFUCKER
    swapV(u, u0); swapV(v, v0);
    advect(2, dt, v, v0, u0, v0, solidCoordinates, solidMap); advect(1, dt, u, u0, u0, v0, solidCoordinates, solidMap);
    project(u, v, p, v0, solidCoordinates, solidMap);
}
void navier::fluid::sStep()
{
    swapV(s, s0); 
    diffuse(0, dt, diff, s, s0, solidCoordinates, solidMap);
    swapV(s, s0); 
    advect(0, dt, s, s0, u, v, solidCoordinates, solidMap);
}
void navier::fluid::step(float dt_)
{
    dt = dt_;
    vStep();
    sStep();
    t += dt;

}

navier::simulation::simulation (fluid Fluid_) : Fluid(Fluid_) 
{
    Nj = Fluid.Nj, Ni = Fluid.Ni, WIDTH = Nj-2, HEIGHT = Ni-2;
}
void navier::simulation::start()
{
    sf::Texture texture;
    sf::Sprite sprite;
    sf::Image image;
    sf::RenderWindow window;
    window.create(sf::VideoMode(WIDTH, HEIGHT, 32), "Navier-Stokes", sf::Style::Default);
    window.setFramerateLimit(fps);
    texture.create(WIDTH, HEIGHT);
    sprite.setTexture(texture);
    sprite.setPosition(0,0);


    std::vector<sf::Uint8> pixels(HEIGHT*WIDTH*4, 0);
    int i, j, ip, jp;
    float t = 0, factorI, factorJ, x , y, x0, y0;
    sf::Vector2i pos0, pos;
    sf::Event event;
    sf::Clock deltaClock;

    initialConditions(Fluid.s, Fluid.p, Fluid.u, Fluid.v);
    while(window.isOpen())
    {
        sf::Time deltaTime = deltaClock.restart();
        const float dt = deltaTime.asSeconds();
        factorI = (float)HEIGHT/(window.getSize().y-1); factorJ = (float)WIDTH/(window.getSize().x-1);
        periodicP(Fluid.p, t, dt); periodicV(Fluid.v, t, dt); periodicU(Fluid.u, t, dt); periodicS(Fluid.s, t, dt); 
        Fluid.step(dt);

        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
            {
                window.close();
            }
            if (event.type == sf::Event::MouseMoved)
            {
                factorI = (float)HEIGHT/(window.getSize().y-1); factorJ = (float)WIDTH/(window.getSize().x-1);
                
            }
            if (sf::Mouse::isButtonPressed(sf::Mouse::Right)) 
            {
                pos = sf::Mouse::getPosition(window);
                y = (pos.y * factorI), x = (pos.x * factorJ);
                ip = (int)y+1; jp = (int)x+1;
                i = ip > 1 ? (ip < Ni-2 ? ip : Ni-2) : 1;	j = jp > 1 ? (jp < Nj-2 ? jp : Nj-2) : 1;	

                Fluid.addV(i, j, radius, mouseV*(x-x0), -mouseV*(y-y0));
                
                pos0 = sf::Mouse::getPosition(window);
                y0 = (pos.y * factorI), x0 = (pos.x * factorJ);
            }
            
        }
    
        if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) 
        {
            pos = sf::Mouse::getPosition(window);
            ip = (int)(pos.y * factorI)+1; jp = (int)(pos.x * factorJ)+1;
            i = ip > 1 ? (ip < Ni-2 ? ip : Ni-2) : 1;	j = jp > 1 ? (jp < Nj-2 ? jp : Nj-2) : 1;	
            Fluid.addS(i, j, radius, dt * mouseS);
        } 
        if (sf::Mouse::isButtonPressed(sf::Mouse::Middle)) 
        {
            pos = sf::Mouse::getPosition(window);
            ip = (int)(pos.y * factorI)+1; jp = (int)(pos.x * factorJ)+1;
            i = ip > 1 ? (ip < Ni-2 ? ip : Ni-2) : 1;	j = jp > 1 ? (jp < Nj-2 ? jp : Nj-2) : 1;	
            Fluid.addP(i, j, radius, dt * mouseP);
        } 

        window.clear(sf::Color::Black);
        switch(draw)
        {
            case 0: paintPixels(Fluid.s, Fluid.solidMap, pixels, r, g, b, a); break;
            case 1: paintPixels(Fluid.u, Fluid.solidMap, pixels, r, g, b, a); break;
            case 2: paintPixels(Fluid.v, Fluid.solidMap, pixels, r, g, b, a); break;
            case 3: paintPixels(Fluid.p, Fluid.solidMap, pixels, r, g, b, a); break;
        }
        texture.update(pixels.data());
        window.draw(sprite);
        window.display();
        t += dt;
    }
}




void navier::advect(int b, float dt, grid &s, grid &s0, grid &u, grid &v, coordinates &solidCoordinates, grid &solidMap)
{
    int i, j, iT, jL, iB, jR, Nj = s.getNj(), Ni = s.getNi(); // assume velocity in cells/s
    float ip, jp, interpX, interpY;
    for (i=1; i<Ni-1; i++)
    {
        for (j=1; j<Nj-1; j++)
        {
            jp = j - dt * u(i, j);  ip = i + dt * v(i, j);
            if (ip < 0.5) ip = 0.5; else if (ip > Ni-1.5) ip = Ni-1.5;
            if (jp < 0.5) jp = 0.5; else if (jp > Nj-1.5) jp = Nj-1.5;
            iT = (int)ip; jL = (int)jp; iB = iT+1; jR = jL+1; interpY = ip - iT; interpX = jp - jL;
            s(i, j) = 
            (s0(iT, jL) * (1 - interpY) + s0(iB, jL) * interpY) * (1 - interpX) + 
            (s0(iT, jR) * (1 - interpY) + s0(iB, jR) * interpY) *      interpX;
        }
    }
    boundaries(b, s);
    lineBoundaries(b, s, solidCoordinates, solidMap);
}
void navier::diffuse(int b, float dt, float diff, grid &s, grid &s0, coordinates &solidCoordinates, grid &solidMap)
{
    int i, j, k, Nj = s.getNj(), Ni = s.getNi();
    float a = dt * diff;

    for (k=0; k<20; k++)
    {
        for (i=1; i<Ni-1; i++)
        {
            for (j=1; j<Nj-1; j++) 
            {
               s(i, j) = (s0(i, j) + a * (s(i+1, j) + s(i-1, j) + s(i, j+1) + s(i, j-1))) / (1 + 4*a);

            }
        }
        boundaries(b, s);
        lineBoundaries(b, s, solidCoordinates, solidMap);
    }

}
void navier::boundaries(int b, grid &s) // 1 for -nx
{ 
    int i, j, Nj = s.getNj(), Ni = s.getNi();
    for (i=1; i < Ni-1; i++)
    {
        s(i, 0) = b==1 ? -s(i, 1): s(i, 1);
        s(i, Nj-1) = b==1 ? -s(i, Nj-2) : s(i, Nj-2);
    }
    for (j=1; j < Nj-1; j++)
    {
        s(0, j) = b==2 ? -s(1, j) : s(1, j);
        s(Ni-1, j) = b==2 ? -s(Ni-2, j) : s(Ni-2, j);
    }
    s(0, 0) = 0.5 * (s(1, 0) + s(0, 1));
    s(Ni-1, Nj-1) = 0.5 * (s(Ni-2, Nj-1) + s(Ni-1, Nj-2));
    s(0, Nj-1) = 0.5 * (s(1, Nj-1) + s(0, Nj-2));
    s(Ni-1, 0) = 0.5 * (s(Ni-1, 1) + s(Ni-2, 0));
}
void navier::lineBoundaries(int b, grid &s, coordinates &solidCoordinates, grid &solidMap)
{
    int N = solidCoordinates.is.size(), i, j, nF;
    float val;
    for (int k=0; k<N; k++)
    {
        i = solidCoordinates.is[k]; j = solidCoordinates.js[k];
        val = 0; nF = 0;
        if (solidMap(i, j+1) != 0) {val += b == 1 ? -s(i, j+1) : s(i, j+1); nF ++;}
        if (solidMap(i, j-1) != 0) {val += b == 1 ? -s(i, j-1) : s(i, j-1); nF ++;}
        if (solidMap(i+1, j) != 0) {val += b == 2 ? -s(i+1, j) : s(i+1, j); nF ++;}
        if (solidMap(i-1, j) != 0) {val += b == 2 ? -s(i-1, j) : s(i-1, j); nF ++;}
        if (nF == 0) //only one fluid corner then
        {
            if (solidMap(i-1, j-1) != 0) {val += 0.5 * (s(i-1, j)+s(i,j-1)) ; nF ++;}
            if (solidMap(i+1, j+1) != 0) {val += 0.5 * (s(i+1, j)+s(i,j+1)) ; nF ++;}
            if (solidMap(i-1, j+1) != 0) {val += 0.5 * (s(i-1, j)+s(i,j+1)) ; nF ++;}
            if (solidMap(i+1, j-1) != 0) {val += 0.5 * (s(i+1, j)+s(i,j-1)) ; nF ++;}       
        }    
        if (nF != 0) {val/=nF;}
        s(i, j) = val;
    }

}
void navier::project(grid &u, grid &v, grid &p, grid &div, coordinates &solidCoordinates, grid &solidMap)
{
    int i, j, k, Nj = u.getNj(), Ni = u.getNi();
    for (i=1; i<Ni-1; i++)
    {
        for (j=1; j<Nj-1; j++)
        {
            div(i, j) = 0.5 * ( u(i, j+1) - u(i, j-1) +
                                        v(i-1, j) - v(i+1, j) ); 
            //p(i, j) = 0;  
        }
    }
    boundaries(0, p); lineBoundaries(0, p, solidCoordinates, solidMap); 
    boundaries(0, div); lineBoundaries(0, div, solidCoordinates, solidMap);
    for (k=0; k<20; k++)
    {
        for (i=1; i<Ni-1; i++)
        {
            for (j=1; j<Nj-1; j++) // p must be initialized
            {
                p(i, j) = (-div(i, j) + p(i+1, j) + p(i-1, j) + p(i, j+1) + p(i, j-1)) * 0.25;

            }
        }
        boundaries(0, p);
        lineBoundaries(0, p, solidCoordinates, solidMap);
    }
    for (i=1; i<Ni-1; i++)
    {
        for (j=1; j<Nj-1; j++) 
        {
            u(i, j) -= 0.5 * (p(i, j+1) - p(i, j-1));
            v(i, j) -= 0.5 * (p(i-1, j) - p(i+1, j));
            
        }
    }
    boundaries(1, u); lineBoundaries(1, u, solidCoordinates, solidMap);
    boundaries(2, v); lineBoundaries(2, v, solidCoordinates, solidMap);
}

void navier::paintPixels(grid &s, grid &solidMap, std::vector<sf::Uint8> &pixels, float  r, float g, float b, int a) // fuck boundaries
{
    int k=0, Nj = s.getNj(), Ni = s.getNi();
    float val;

    for(int i=1; i<Ni-1; i++)
    {
        for (int j=1; j<Nj-1; j++)
        {
            if (solidMap(i, j) != 0)
            {
                pixels[k  ] = 255;
                pixels[k+1] = 255;
                pixels[k+2] = 255;
                pixels[k+3] = 255; 
            }
            else
            {
                val = s(i, j);
                pixels[k  ] = val > 255 ? 255 : (sf::Uint8)val * r;
                pixels[k+1] = val > 255 ? 255 : (sf::Uint8)val * g;
                pixels[k+2] = val > 255 ? 255 : (sf::Uint8)val * b;
                pixels[k+3] =                                    a;
            }
            k += 4;
        }
    }
}
void navier::addSource(int radius, grid &s, int i, int j, int b, float val, coordinates &solidCoordinates, grid &solidMap) // assume solid walls, Ni is the whole Ni, i and j too 
{
    int Nj = s.getNj(), Ni = s.getNi();
    int spaceBot = (Ni-2) - i;
    int spaceTop = i - 1;
    int spaceRight = (Nj-2) - j;
    int spaceLeft = j - 1;
    int radBot   = spaceBot   > radius ? radius :   spaceBot;
    int radTop   = spaceTop   > radius ? radius :   spaceTop;
    int radLeft  = spaceLeft  > radius ? radius :  spaceLeft;
    int radRight = spaceRight > radius ? radius : spaceRight;
    for (int ik = i-radTop; ik <= i + radBot; ik++)
    {
        for (int jk = j-radLeft; jk <= j + radRight; jk++)
        {
            s(ik, jk) += val;
            //std::cout << "val is now "<< S(ik, jk)  << "\n";
        } 
    }
    boundaries(b, s);
    lineBoundaries(b, s, solidCoordinates, solidMap);
}