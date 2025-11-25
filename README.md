# Navier-Stokes Fluid Simulation - Refactored

A professional C++ implementation of incompressible fluid dynamics using the Navier-Stokes equations with support for solid obstacles and real-time interaction.

## Architecture Overview

The project has been refactored following modern C++ best practices and SOLID principles:

### Core Components

1. **Grid** (`Grid.h/cpp`)
   - 2D container for fluid properties
   - Efficient row-major storage
   - Bounds checking and utility methods
   - Drawing primitives for obstacles

2. **FluidSolver** (`FluidSolver.h/cpp`)
   - Core Navier-Stokes solver
   - Implements operator splitting method:
     - Advection (semi-Lagrangian)
     - Diffusion (implicit Gauss-Seidel)
     - Projection (Helmholtz decomposition)
   - Solid boundary handling
   - Vorticity computation

3. **Renderer** (`Renderer.h/cpp`)
   - Visualization layer (SFML-based)
   - Multiple color mapping modes
   - Gradient and solid color visualization
   - Screen-grid coordinate conversion

4. **Simulation** (`Simulation.h/cpp`)
   - Main orchestrator
   - Event handling and user interaction
   - Mouse-based fluid manipulation
   - Keyboard shortcuts for visualization modes

5. **Config** (`Config.h`)
   - Centralized configuration constants
   - Physical parameters
   - Visualization settings
   - Interaction defaults

## Key Improvements

### Design Patterns
✅ **Separation of Concerns**: Solver, rendering, and interaction are independent  
✅ **Single Responsibility**: Each class has one clear purpose  
✅ **Dependency Inversion**: Using std::function instead of raw function pointers  
✅ **Encapsulation**: Private members with proper accessors  
✅ **RAII**: Automatic resource management with smart pointers  

### Code Quality
✅ **Modern C++**: std::function, std::clamp, range-based loops, structured bindings  
✅ **Const Correctness**: Proper const methods and parameters  
✅ **Clear Naming**: Descriptive variable and method names (no more single letters)  
✅ **Documentation**: Doxygen-style comments  
✅ **Type Safety**: enum class instead of integers  
✅ **Error Handling**: Bounds checking and exception handling  

### Performance
✅ **Memory Efficiency**: Contiguous storage, cache-friendly access  
✅ **Reduced Copying**: Using references and move semantics  
✅ **Optimized Loops**: Iterator-based where appropriate  

## Building the Project

### Requirements
- C++17 or later
- SFML 2.5+
- CMake (optional, for build system)

### Compilation
```bash
g++ -std=c++17 -o fluid_sim \
    main.cpp \
    Grid.cpp \
    FluidSolver.cpp \
    Renderer.cpp \
    Simulation.cpp \
    -lsfml-graphics -lsfml-window -lsfml-system
```

## Usage

### Basic Example
```cpp
#include "Simulation.h"

int main() {
    // Create 100x100 grid with 4x screen scaling
    navier::Simulation sim(100, 100, 4.0f);
    
    // Configure physics
    sim.setPhysicalParameters(1.0f, 1.0f, 0.005f);
    
    // Add obstacles
    sim.addObstacles({
        navier::Line(40, 40, 20, true, 2)
    });
    
    // Run simulation
    sim.run();
    
    return 0;
}
```

### Interactive Controls
- **Mouse**: Click and drag to add scalar/velocity/pressure
- **0**: View scalar field (dye)
- **1**: View horizontal velocity (u)
- **2**: View vertical velocity (v)
- **3**: View pressure field
- **4**: View vorticity (curl)
- **+/-**: Increase/decrease brush radius
- **G**: Toggle gradient visualization
- **N**: Toggle value normalization

## Advanced Features

### Periodic Forcing
```cpp
sim.setPeriodicForcing(
    [](Grid& u, double t, double dt) {
        // Add time-varying forces
        u(50, 50) += std::sin(t) * 100.0f;
    },
    nullptr, nullptr, nullptr
);
```

### Custom Boundaries
```cpp
sim.setBoundaryConditions(
    BoundaryType::HORIZONTAL,  // Reflect horizontal velocity at walls
    BoundaryType::VERTICAL,    // Reflect vertical velocity at walls
    BoundaryType::NONE,        // Free-slip for scalar
    BoundaryType::NONE         // Free-slip for pressure
);
```

### Visualization Configuration
```cpp
ColorMap map;
map.r = 1.0f; map.g = 0.0f; map.b = 0.0f;
map.useGradient = false;
map.normalize = true;
sim.setColorMap(map);
```

## Physical Model

The solver implements the incompressible Navier-Stokes equations:

```
∂u/∂t + (u·∇)u = -∇p + ν∇²u + f
∇·u = 0
```

Where:
- **u**: velocity field
- **p**: pressure field
- **ν**: kinematic viscosity
- **f**: external forces

### Numerical Method
1. **Advection**: Semi-Lagrangian backtrace with bilinear interpolation
2. **Diffusion**: Implicit solver using Gauss-Seidel iteration
3. **Projection**: Helmholtz decomposition to enforce incompressibility
4. **Solid Boundaries**: No-slip/free-slip conditions at obstacles

## Project Structure

```
navier-stokes-main/
├── Config.h              # Configuration constants
├── Grid.h/cpp            # 2D grid container
├── FluidSolver.h/cpp     # Navier-Stokes solver
├── Renderer.h/cpp        # Visualization layer
├── Simulation.h/cpp      # Main simulation controller
├── main.cpp              # Example application
└── README.md             # This file
```

## Performance Considerations

- **Grid Size**: 100x100 typical, can handle up to 500x500 on modern hardware
- **Time Step**: Automatically clamped for stability (0.001-0.1)
- **Iterations**: 20 Gauss-Seidel iterations per solve (configurable)
- **FPS**: Default 60, adjustable based on performance needs

## Future Enhancements

Potential improvements for the future:
- [ ] Multigrid solver for faster convergence
- [ ] GPU acceleration (CUDA/OpenCL)
- [ ] Adaptive time stepping (CFL condition)
- [ ] Save/load simulation states
- [ ] More obstacle shapes (circles, arbitrary polygons)
- [ ] Turbulence modeling
- [ ] Multi-threaded solver

## References

- **Jos Stam**: "Real-Time Fluid Dynamics for Games" (2003)
- **Bridson**: "Fluid Simulation for Computer Graphics" (2015)
- **Griebel et al.**: "Numerical Simulation in Fluid Dynamics" (1998)

## License

MIT License

## Credits

Original implementation refactored with modern C++ practices and professional architecture patterns.
