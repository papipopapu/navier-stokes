#ifndef SIMULATION_H
#define SIMULATION_H

#include "FluidSolver.h"
#include "Renderer.h"
#include "Config.h"
#include <SFML/Graphics.hpp>

namespace navier {

/**
 * @brief Main simulation controller
 * Orchestrates the fluid solver, renderer, and user interaction
 */
class Simulation {
private:
    FluidSolver m_solver;
    Renderer m_renderer;
    
    // Simulation parameters
    float m_playSpeed;
    float m_timeStep;
    uint8_t m_targetFPS;
    
    // Interaction parameters
    int m_brushRadius;
    float m_scalarStrength;
    float m_velocityStrength;
    float m_pressureStrength;
    
    // Maximum addition values
    float m_maxAddScalar;
    float m_maxAddVelocity;
    float m_maxAddPressure;
    
    // Mouse state
    sf::Vector2i m_previousMousePos;
    bool m_isMousePressed;
    
    // Periodic forcing callbacks
    FluidSolver::PeriodicFunction m_periodicU;
    FluidSolver::PeriodicFunction m_periodicV;
    FluidSolver::PeriodicFunction m_periodicScalar;
    FluidSolver::PeriodicFunction m_periodicPressure;

public:
    Simulation(size_t rows, size_t cols, float screenScale = 1.0f);
    
    // Main simulation loop
    void run();
    
    // Solver configuration
    void setPhysicalParameters(float viscosity, float diffusion, float dissipation);
    void setBoundaryConditions(BoundaryType u, BoundaryType v, 
                              BoundaryType scalar, BoundaryType pressure);
    void addObstacles(std::initializer_list<Line> lines);
    
    // Visualization configuration
    void setColorMap(const ColorMap& colorMap);
    void setVisualizationMode(VisualizationConfig::DrawMode mode);
    void setDrawRange(float minVal, float maxVal);
    
    // Interaction configuration
    void setInteractionStrengths(float scalar, float velocity, float pressure);
    void setMaxAdditions(float scalar, float velocity, float pressure);
    void setBrushRadius(int radius);
    
    // Time stepping
    void setPlaySpeed(float speed) { m_playSpeed = speed; }
    void setTargetFPS(uint8_t fps);
    
    // Periodic forcing
    void setPeriodicForcing(FluidSolver::PeriodicFunction u,
                           FluidSolver::PeriodicFunction v,
                           FluidSolver::PeriodicFunction scalar,
                           FluidSolver::PeriodicFunction pressure);
    
    // Initial conditions
    void setInitialConditions(FluidSolver::InitialConditionFunction initFunc);
    
    // Access to solver (for advanced use)
    FluidSolver& solver() { return m_solver; }
    const FluidSolver& solver() const { return m_solver; }

private:
    void handleEvents();
    void handleMouseInput();
    void update(float dt);
    void render();
    
    void processKeyPress(sf::Keyboard::Key key);
    void addSourceAtMouse(int gridI, int gridJ);
    
    const Grid& getCurrentVisualizationField() const;
};

} // namespace navier

#endif // SIMULATION_H
