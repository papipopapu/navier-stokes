![Evolution of a simulated ring vortex](images/timeline.jpg)

# Navier-Stokes Fluid Simulation 

A robust C++ implementation of incompressible fluid dynamics using the Navier-Stokes equations with support for solid obstacles and real-time interaction. Created to study the behavior of ring vortices.

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


## Building the Project

### Requirements
- C++17 or later
- SFML 2.5+
- CMake (optional, for build system)

### Compilation
```bash
g++ -std=c++17 -o fluid_sim \
    core/main.cpp \
    core/Grid.cpp \
    core/FluidSolver.cpp \
    core/Renderer.cpp \
    core/Simulation.cpp \
    -I core \
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

## Solver Theory

### Motivation and Physical Background

The study of vortex dynamics and structure in fluids is crucial for weather prediction and engineering design of systems operating in fluid environments. Given the complexity of reality, theoretical models based on empirical methods are necessary to study turbulent flow behavior and find regular patterns and predictable trends in nature.

Von Kármán vortex streets and ring vortices are regularly observed in natural phenomena, from atmospheric patterns around the Canary Islands to historical events. This project aims to explain such phenomena through simplified 2D simulations.

A vortex is a region of rotating fluid flow. According to Bernoulli's principle, the pressure in a vortex is lower in the central region where the velocity is higher, increasing as we move away from the center.

Ring vortices arise when fluid is propelled through an opening or impulsively, forming characteristic toroidal (doughnut-shaped) structures. The vorticity at the boundary of the jet or puff rolls up into a ring due to velocity gradients, creating a self-propagating vortex with circulation around the ring's cross-section. Von Kármán vortex streets are repetitive patterns of swirling vortices caused by vortex shedding, responsible for unstable boundary layer separation in fluid flow around solids.

### Computational Approach

The simulation solves the incompressible Navier-Stokes equations for both the velocity field:

```
∂u/∂t = -(u·∇)u - (1/ρ)∇p + ν∇²u + f     (Eq. 1)
```

(In the non-dimensional form or when ρ=1, this simplifies to: ∂u/∂t + (u·∇)u = -∇p + ν∇²u + f)

and for scalar fields (e.g., dye density):

```
∂s/∂t = -(u·∇)s + κ∇²s + S                (Eq. 2)
```

The implementation uses an **Eulerian** approach, computing relevant quantities at fixed points rather than tracking particle trajectories (Lagrangian Particle Tracking). All methods are based on finite differences, working with a discretized environment that provides extensive equation systems.

### MAC Grid Structure

The simulation employs a **MAC (Marker and Cell)** grid structure:
- Scalar quantities are computed and stored at cell centers
- Vector components are stored at cell faces (u-component on left/right faces, v-component on top/bottom faces)
- This staggered arrangement improves accuracy of discretized operators
- Uniform distribution simplifies operator effects and enables value interpolation at any point in the fluid

### Operator Splitting Method

Each iteration can be viewed as the composition of operators:

```
u^(n+1) = P ∘ F ∘ D ∘ A(u^n)       (velocity)
s^(n+1) = S ∘ D ∘ A(s^n)           (scalar)
```

Where:
- **A (Advection)**: Represents the term (u·∇), the fluid moving through itself
- **D (Diffusion)**: Corresponds to the Laplacian term ∇²u
- **F (Forces)** and **S (Sources)**: Addition of external forces and scalar sources
- **P (Projection)**: Applies pressure gradient and enforces incompressibility

#### Advection (Operator A)

Uses a **semi-Lagrangian backward trace** for stability. For each cell and each physical quantity:

```
X_prev = X_cell - u^old · Δt
Q^new = interpolate(X_prev)
```

This backward tracing is more stable than forward methods. The implementation supports both Euler and Runge-Kutta 4 methods.

#### Diffusion (Operator D)

To maintain stability, an **implicit formulation** is used:

```
Q^old = Q^new + Δt · ν∇²Q^new
```

where ν is the diffusion coefficient (kinematic viscosity for velocity, or diffusivity κ for scalar fields). This creates a system with as many equations as fluid cells, solved using the **Gauss-Seidel relaxation method**. Note that numerical diffusion often arises naturally and may not always need explicit treatment.

#### Projection (Operator P)

The projection step introduces the pressure gradient term and enforces incompressibility using the **Helmholtz-Hodge decomposition theorem**, which decomposes any vector field into an irrotational and an incompressible component:

```
w = ∇Φ + curl(A)
```

Considering that after applying all operators except P, the velocity field may be compressible, we associate the desired velocity with the incompressible component:

```
w ≡ u_old = ∇p + u_new
u_new = u_old - ∇p  ⟹  ∇²p = ∇·u_old   (Eq. 3)
```

Solving this **Poisson equation** (Eq. 3) yields the pressure field p that makes the velocity incompressible. The system is solved both by Gauss-Seidel iteration and by direct sparse matrix inversion (more accurate but slower).

#### Boundary Conditions

Throughout this process, boundary conditions are applied at cells bordering solids or simulation boundaries. For solid walls, the perpendicular velocity component at the boundary is set to zero.

## Ring Vortex Results

Contrasting with basic theory about ring vortices, which provides the expression:

```
U = (πω₀a²)/(4πR) · (ln(8R/a) - 1/4)
```

where:
- **a**: characteristic vortex core radius
- **R**: toroidal radius
- **ω₀**: vorticity

Our simulations deviate from this theoretical regime. The model assumes uniform vorticity in each of the two small vortices and that their radius is much smaller than the toroidal radius R (a << R). In our 2D simulations, these approximations are not reasonable.

### Key Findings

The most studied quantity is **propagation velocity**, which shows close relationships with other quantities:

1. **Velocity vs. Energy**: Approximately **linear relationship** between translation velocity and initial energy input

2. **Velocity vs. Mean Radius**: Approximately **linear relationship** between translation velocity and mean vortex radius. In our observations, the vortex cross-section has inner and outer radii forming a ring, and we measure the mean radius (a+b)/2 where a and b are the inner and outer radii respectively.

   These linear relationships indicate that velocity, energy, and radius are closely interconnected. If we assume the characteristic radius a and toroidal radius R are proportional (R = k·a), the theoretical formula becomes U = (πω₀a)/(4k)·(ln(8k) - 1/4), which is indeed linear with a, though this remains conjectural.

3. **Velocity vs. Vorticity**: The relationship shows approximately **U ∝ √ω**, where ω is vorticity and U is translation velocity

### Conclusion

While results do not conform to the simplified theoretical model, this is expected given that our system does not meet the approximations assumed by that model. The relationships obtained escape this simplification, showing behavior more complex than the idealized case.

## Project Structure

```
navier-stokes-main/
├── core/
│   ├── Config.h          # Configuration constants
│   ├── Grid.h/cpp        # 2D grid container
│   ├── FluidSolver.h/cpp # Navier-Stokes solver
│   ├── Renderer.h/cpp    # Visualization layer
│   ├── Simulation.h/cpp  # Main simulation controller
│   └── main.cpp          # Example application
├── images/
│   └── timeline.jpg      # Simulation timeline image
├── CMakeLists.txt        # Build configuration
└── README.md             # This file
```

## Performance Considerations

- **Grid Size**: 100x100 typical, can handle up to 500x500 on modern hardware
- **Time Step**: Automatically clamped for stability (0.001-0.1)
- **Iterations**: 20 Gauss-Seidel iterations per solve (configurable)
- **FPS**: Default 60, adjustable based on performance needs

## License

MIT License






