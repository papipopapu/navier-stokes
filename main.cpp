#include "Simulation.h"
#include <iostream>

using namespace navier;

int main() {
    try {
        // Create simulation with 100x100 grid, 4x screen scale
        Simulation sim(100, 100, 4.0f);
        
        // Configure physical parameters
        sim.setPhysicalParameters(
            1.0f,    // viscosity
            1.0f,    // diffusion
            0.005f   // dissipation
        );
        
        // Set boundary conditions
        sim.setBoundaryConditions(
            BoundaryType::HORIZONTAL,  // u velocity
            BoundaryType::VERTICAL,    // v velocity
            BoundaryType::NONE,        // scalar
            BoundaryType::NONE         // pressure
        );
        
        // Add obstacles (example: box in the center)
        sim.addObstacles({
            Line(40, 40, 20, true, 2),   // Left wall
            Line(40, 58, 20, true, 2),   // Right wall
            Line(40, 40, 20, false, 2),  // Top wall
            Line(58, 40, 20, false, 2)   // Bottom wall
        });
        
        // Configure visualization
        ColorMap colorMap;
        colorMap.useGradient = true;
        colorMap.normalize = true;
        sim.setColorMap(colorMap);
        
        // Configure interaction
        sim.setBrushRadius(3);
        sim.setInteractionStrengths(2000.0f, 10.0f, 100000.0f);
        
        // Set simulation speed and FPS
        sim.setPlaySpeed(0.2f);
        sim.setTargetFPS(60);
        
        // Print controls
        std::cout << "=== Navier-Stokes Fluid Simulation ===\n";
        std::cout << "Controls:\n";
        std::cout << "  Mouse: Click and drag to add forces/scalar\n";
        std::cout << "  0: View scalar field\n";
        std::cout << "  1: View velocity U\n";
        std::cout << "  2: View velocity V\n";
        std::cout << "  3: View pressure\n";
        std::cout << "  4: View vorticity\n";
        std::cout << "  +/-: Increase/decrease brush radius\n";
        std::cout << "  G: Toggle gradient mode\n";
        std::cout << "  N: Toggle normalization\n";
        std::cout << "  ESC: Exit\n\n";
        
        // Run the simulation
        sim.run();
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}