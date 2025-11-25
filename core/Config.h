#ifndef CONFIG_H
#define CONFIG_H

#include <cstdint>

namespace navier {

// Simulation Constants
struct SimulationConfig {
    // Physical parameters
    static constexpr float DEFAULT_VISCOSITY = 1.0f;
    static constexpr float DEFAULT_DIFFUSION = 1.0f;
    static constexpr float DEFAULT_DISSIPATION = 0.005f;
    
    // Solver parameters
    static constexpr int GAUSS_SEIDEL_ITERATIONS = 20;
    static constexpr float RELAXATION_FACTOR = 1.9f;
    static constexpr float CONVERGENCE_THRESHOLD = 1e-6f;
    
    // Time stepping
    static constexpr float MIN_TIMESTEP = 0.001f;
    static constexpr float MAX_TIMESTEP = 0.1f;
    
    // Boundary conditions
    static constexpr int BOUNDARY_NONE = 0;
    static constexpr int BOUNDARY_HORIZONTAL = 1;
    static constexpr int BOUNDARY_VERTICAL = 2;
};

// Visualization Constants
struct VisualizationConfig {
    // Default colors
    static constexpr float DEFAULT_R = 0.0f;
    static constexpr float DEFAULT_G = 0.0f;
    static constexpr float DEFAULT_B = 1.0f;
    static constexpr uint8_t DEFAULT_ALPHA = 255;
    
    // Visualization modes
    enum class DrawMode : uint8_t {
        SCALAR = 0,
        VELOCITY_U = 1,
        VELOCITY_V = 2,
        PRESSURE = 3,
        VORTICITY = 4
    };
    
    // FPS limits
    static constexpr uint8_t DEFAULT_FPS = 60;
    static constexpr uint8_t MIN_FPS = 10;
    static constexpr uint8_t MAX_FPS = 120;
    
    // Typical value ranges for normalization
    static constexpr float MAX_VELOCITY = 400.0f;
    static constexpr float MAX_SCALAR = 200.0f;
    static constexpr float MAX_PRESSURE = 1000.0f;
};

// Interaction Constants
struct InteractionConfig {
    // Mouse interaction
    static constexpr int DEFAULT_BRUSH_RADIUS = 3;
    static constexpr int MIN_BRUSH_RADIUS = 1;
    static constexpr int MAX_BRUSH_RADIUS = 20;
    
    // Default interaction strengths
    static constexpr float DEFAULT_SCALAR_STRENGTH = 2000.0f;
    static constexpr float DEFAULT_VELOCITY_STRENGTH = 10.0f;
    static constexpr float DEFAULT_PRESSURE_STRENGTH = 100000.0f;
    
    // Maximum addition values
    static constexpr float MAX_ADD_SCALAR = 2000.0f;
    static constexpr float MAX_ADD_VELOCITY = 100.0f;
    static constexpr float MAX_ADD_PRESSURE = 100000.0f;
};

} // namespace navier

#endif // CONFIG_H
