#include "Simulation.h"
#include <iostream>

namespace navier {

Simulation::Simulation(size_t rows, size_t cols, float screenScale)
    : m_solver(rows, cols)
    , m_renderer(rows, cols, screenScale)
    , m_playSpeed(0.2f)
    , m_timeStep(0.016f) // ~60 FPS default
    , m_targetFPS(VisualizationConfig::DEFAULT_FPS)
    , m_brushRadius(InteractionConfig::DEFAULT_BRUSH_RADIUS)
    , m_scalarStrength(InteractionConfig::DEFAULT_SCALAR_STRENGTH)
    , m_velocityStrength(InteractionConfig::DEFAULT_VELOCITY_STRENGTH)
    , m_pressureStrength(InteractionConfig::DEFAULT_PRESSURE_STRENGTH)
    , m_maxAddScalar(InteractionConfig::MAX_ADD_SCALAR)
    , m_maxAddVelocity(InteractionConfig::MAX_ADD_VELOCITY)
    , m_maxAddPressure(InteractionConfig::MAX_ADD_PRESSURE)
    , m_isMousePressed(false)
{
    m_renderer.window().setFramerateLimit(m_targetFPS);
}

void Simulation::run() {
    sf::Clock clock;
    
    while (m_renderer.isOpen()) {
        handleEvents();
        
        // Calculate frame time
        float frameTime = clock.restart().asSeconds();
        float simulationTime = frameTime * m_playSpeed;
        
        // Update simulation
        update(simulationTime);
        
        // Render
        render();
    }
}

void Simulation::handleEvents() {
    sf::Event event;
    while (m_renderer.window().pollEvent(event)) {
        if (event.type == sf::Event::Closed) {
            m_renderer.close();
        }
        else if (event.type == sf::Event::KeyPressed) {
            processKeyPress(event.key.code);
        }
        else if (event.type == sf::Event::MouseButtonPressed) {
            if (event.mouseButton.button == sf::Mouse::Left) {
                m_isMousePressed = true;
                m_previousMousePos = sf::Vector2i(event.mouseButton.x, event.mouseButton.y);
            }
        }
        else if (event.type == sf::Event::MouseButtonReleased) {
            if (event.mouseButton.button == sf::Mouse::Left) {
                m_isMousePressed = false;
            }
        }
    }
    
    handleMouseInput();
}

void Simulation::handleMouseInput() {
    if (m_isMousePressed) {
        sf::Vector2i mousePos = m_renderer.getMousePosition();
        int gridI, gridJ;
        m_renderer.screenToGrid(mousePos.x, mousePos.y, gridI, gridJ);
        
        // Check bounds
        if (gridI >= 0 && gridI < static_cast<int>(m_solver.rows()) &&
            gridJ >= 0 && gridJ < static_cast<int>(m_solver.cols())) {
            addSourceAtMouse(gridI, gridJ);
        }
        
        m_previousMousePos = mousePos;
    }
}

void Simulation::addSourceAtMouse(int gridI, int gridJ) {
    auto mode = m_renderer.getDrawMode();
    
    switch (mode) {
        case VisualizationConfig::DrawMode::SCALAR:
            m_solver.addScalar(gridI, gridJ, m_brushRadius, m_scalarStrength, m_maxAddScalar);
            break;
            
        case VisualizationConfig::DrawMode::VELOCITY_U:
        case VisualizationConfig::DrawMode::VELOCITY_V: {
            // Calculate velocity from mouse movement
            sf::Vector2i currentPos = m_renderer.getMousePosition();
            sf::Vector2i delta = currentPos - m_previousMousePos;
            
            float u = static_cast<float>(delta.x) * m_velocityStrength;
            float v = static_cast<float>(delta.y) * m_velocityStrength;
            
            m_solver.addVelocity(gridI, gridJ, m_brushRadius, u, v, 
                                m_maxAddVelocity, m_maxAddVelocity);
            break;
        }
        
        case VisualizationConfig::DrawMode::PRESSURE:
            m_solver.addPressure(gridI, gridJ, m_brushRadius, m_pressureStrength, m_maxAddPressure);
            break;
            
        default:
            break;
    }
}

void Simulation::update(float dt) {
    // Clamp timestep to reasonable range
    dt = std::clamp(dt, SimulationConfig::MIN_TIMESTEP, SimulationConfig::MAX_TIMESTEP);
    
    // Apply periodic forcing if set
    if (m_periodicU || m_periodicV || m_periodicScalar || m_periodicPressure) {
        m_solver.applyPeriodicForcing(m_periodicU, m_periodicV, 
                                     m_periodicScalar, m_periodicPressure, dt);
    }
    
    // Step the simulation
    m_solver.step(dt);
}

void Simulation::render() {
    const Grid& field = getCurrentVisualizationField();
    m_renderer.render(field, m_solver.solidMap());
    m_renderer.display();
}

const Grid& Simulation::getCurrentVisualizationField() const {
    switch (m_renderer.getDrawMode()) {
        case VisualizationConfig::DrawMode::VELOCITY_U:
            return m_solver.velocityU();
        case VisualizationConfig::DrawMode::VELOCITY_V:
            return m_solver.velocityV();
        case VisualizationConfig::DrawMode::PRESSURE:
            return m_solver.pressure();
        case VisualizationConfig::DrawMode::VORTICITY:
            return m_solver.vorticity();
        case VisualizationConfig::DrawMode::SCALAR:
        default:
            return m_solver.scalar();
    }
}

void Simulation::processKeyPress(sf::Keyboard::Key key) {
    switch (key) {
        case sf::Keyboard::Num0:
            m_renderer.setDrawMode(VisualizationConfig::DrawMode::SCALAR);
            m_renderer.setDrawRange(0.0f, VisualizationConfig::MAX_SCALAR);
            std::cout << "Visualization: Scalar field\n";
            break;
            
        case sf::Keyboard::Num1:
            m_renderer.setDrawMode(VisualizationConfig::DrawMode::VELOCITY_U);
            m_renderer.setDrawRange(-VisualizationConfig::MAX_VELOCITY, 
                                   VisualizationConfig::MAX_VELOCITY);
            std::cout << "Visualization: Velocity U\n";
            break;
            
        case sf::Keyboard::Num2:
            m_renderer.setDrawMode(VisualizationConfig::DrawMode::VELOCITY_V);
            m_renderer.setDrawRange(-VisualizationConfig::MAX_VELOCITY, 
                                   VisualizationConfig::MAX_VELOCITY);
            std::cout << "Visualization: Velocity V\n";
            break;
            
        case sf::Keyboard::Num3:
            m_renderer.setDrawMode(VisualizationConfig::DrawMode::PRESSURE);
            m_renderer.setDrawRange(-VisualizationConfig::MAX_PRESSURE, 
                                   VisualizationConfig::MAX_PRESSURE);
            std::cout << "Visualization: Pressure\n";
            break;
            
        case sf::Keyboard::Num4:
            m_renderer.setDrawMode(VisualizationConfig::DrawMode::VORTICITY);
            std::cout << "Visualization: Vorticity\n";
            break;
            
        case sf::Keyboard::Add:
        case sf::Keyboard::Equal:
            m_brushRadius = std::min(m_brushRadius + 1, InteractionConfig::MAX_BRUSH_RADIUS);
            std::cout << "Brush radius: " << m_brushRadius << "\n";
            break;
            
        case sf::Keyboard::Subtract:
        case sf::Keyboard::Hyphen:
            m_brushRadius = std::max(m_brushRadius - 1, InteractionConfig::MIN_BRUSH_RADIUS);
            std::cout << "Brush radius: " << m_brushRadius << "\n";
            break;
            
        case sf::Keyboard::G:
            // Toggle gradient mode
            {
                ColorMap map;
                map.useGradient = true;
                m_renderer.setGradientMode(true);
                std::cout << "Gradient mode: ON\n";
            }
            break;
            
        case sf::Keyboard::N:
            // Toggle normalization
            m_renderer.setNormalization(true);
            std::cout << "Normalization: ON\n";
            break;
            
        case sf::Keyboard::Space:
            // Reset simulation
            std::cout << "Reset simulation (not implemented yet)\n";
            break;
            
        default:
            break;
    }
}

// Configuration methods
void Simulation::setPhysicalParameters(float viscosity, float diffusion, float dissipation) {
    m_solver.setViscosity(viscosity);
    m_solver.setDiffusion(diffusion);
    m_solver.setDissipation(dissipation);
}

void Simulation::setBoundaryConditions(BoundaryType u, BoundaryType v, 
                                      BoundaryType scalar, BoundaryType pressure) {
    m_solver.setBoundaryConditions(u, v, scalar, pressure);
}

void Simulation::addObstacles(std::initializer_list<Line> lines) {
    m_solver.addObstacleLines(lines);
}

void Simulation::setColorMap(const ColorMap& colorMap) {
    m_renderer.setColorMap(colorMap);
}

void Simulation::setVisualizationMode(VisualizationConfig::DrawMode mode) {
    m_renderer.setDrawMode(mode);
}

void Simulation::setDrawRange(float minVal, float maxVal) {
    m_renderer.setDrawRange(minVal, maxVal);
}

void Simulation::setInteractionStrengths(float scalar, float velocity, float pressure) {
    m_scalarStrength = scalar;
    m_velocityStrength = velocity;
    m_pressureStrength = pressure;
}

void Simulation::setMaxAdditions(float scalar, float velocity, float pressure) {
    m_maxAddScalar = scalar;
    m_maxAddVelocity = velocity;
    m_maxAddPressure = pressure;
}

void Simulation::setBrushRadius(int radius) {
    m_brushRadius = std::clamp(radius, 
                              InteractionConfig::MIN_BRUSH_RADIUS, 
                              InteractionConfig::MAX_BRUSH_RADIUS);
}

void Simulation::setTargetFPS(uint8_t fps) {
    m_targetFPS = std::clamp(fps, 
                            VisualizationConfig::MIN_FPS, 
                            VisualizationConfig::MAX_FPS);
    m_renderer.window().setFramerateLimit(m_targetFPS);
}

void Simulation::setPeriodicForcing(FluidSolver::PeriodicFunction u,
                                   FluidSolver::PeriodicFunction v,
                                   FluidSolver::PeriodicFunction scalar,
                                   FluidSolver::PeriodicFunction pressure) {
    m_periodicU = u;
    m_periodicV = v;
    m_periodicScalar = scalar;
    m_periodicPressure = pressure;
}

void Simulation::setInitialConditions(FluidSolver::InitialConditionFunction initFunc) {
    // Not implemented in this refactor - would need to expose grids
    // Could be added if needed
}

} // namespace navier
