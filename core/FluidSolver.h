#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

#include "Grid.h"
#include "Config.h"
#include <functional>
#include <initializer_list>

namespace navier {

/**
 * @brief Boundary condition types for fluid properties
 */
enum class BoundaryType {
    NONE = 0,        // No special boundary treatment
    HORIZONTAL = 1,  // Reflect horizontal velocity component
    VERTICAL = 2     // Reflect vertical velocity component
};

/**
 * @brief Core Navier-Stokes fluid solver
 * Implements incompressible fluid simulation using operator splitting
 */
class FluidSolver {
public:
    // Type aliases for callbacks
    using PeriodicFunction = std::function<void(Grid&, double time, double dt)>;
    using InitialConditionFunction = std::function<void(Grid& scalar, Grid& pressure, Grid& u, Grid& v)>;

private:
    // Grid dimensions
    size_t m_rows;
    size_t m_cols;
    
    // Simulation time
    double m_time;
    
    // Velocity fields (u = horizontal, v = vertical)
    Grid m_velocityU;
    Grid m_velocityV;
    Grid m_velocityU_prev;
    Grid m_velocityV_prev;
    
    // Scalar field (e.g., dye, temperature)
    Grid m_scalar;
    Grid m_scalar_prev;
    
    // Pressure field
    Grid m_pressure;
    
    // Vorticity field (curl of velocity)
    Grid m_vorticity;
    
    // Solid obstacle map (0 = fluid, non-zero = solid)
    Grid m_solidMap;
    SolidCoordinates m_solidCoordinates;
    
    // Physical parameters
    float m_viscosity;      // Kinematic viscosity
    float m_diffusion;      // Scalar diffusion coefficient
    float m_dissipation;    // Dissipation rate
    
    // Boundary conditions
    BoundaryType m_boundaryU;
    BoundaryType m_boundaryV;
    BoundaryType m_boundaryScalar;
    BoundaryType m_boundaryPressure;

public:
    FluidSolver(size_t rows, size_t cols);
    
    // Physical parameter setters
    void setViscosity(float viscosity) { m_viscosity = viscosity; }
    void setDiffusion(float diffusion) { m_diffusion = diffusion; }
    void setDissipation(float dissipation) { m_dissipation = dissipation; }
    
    // Boundary condition setters
    void setBoundaryConditions(BoundaryType u, BoundaryType v, 
                              BoundaryType scalar, BoundaryType pressure);
    
    // Obstacle management
    void addObstacleLine(const Line& line);
    void addObstacleLines(std::initializer_list<Line> lines);
    void clearObstacles();
    
    // Source addition (for interaction)
    void addScalar(int i, int j, int radius, float value, float maxValue);
    void addPressure(int i, int j, int radius, float value, float maxValue);
    void addVelocity(int i, int j, int radius, float u, float v, float maxU, float maxV);
    
    // Time stepping
    void step(float dt);
    void velocityStep(float dt);
    void scalarStep(float dt);
    
    // Apply periodic forcing (optional)
    void applyPeriodicForcing(PeriodicFunction forceU, PeriodicFunction forceV,
                            PeriodicFunction forceScalar, PeriodicFunction forcePressure, float dt);
    
    // Accessors (const)
    const Grid& velocityU() const { return m_velocityU; }
    const Grid& velocityV() const { return m_velocityV; }
    const Grid& scalar() const { return m_scalar; }
    const Grid& pressure() const { return m_pressure; }
    const Grid& vorticity() const { return m_vorticity; }
    const Grid& solidMap() const { return m_solidMap; }
    
    double time() const { return m_time; }
    size_t rows() const { return m_rows; }
    size_t cols() const { return m_cols; }
    
    // Compute derived quantities
    void computeVorticity();

private:
    // Core solver methods
    void advect(BoundaryType boundary, float dt, Grid& field, const Grid& field_prev);
    void diffuse(BoundaryType boundary, float dt, float diffCoeff, Grid& field, const Grid& field_prev);
    void project();
    void applyDissipation(Grid& field, float dt);
    
    // Boundary handling
    void applyBoundaryConditions(BoundaryType boundary, Grid& field);
    void applySolidBoundaries(BoundaryType boundary, Grid& field);
    
    // Linear solver for diffusion and pressure projection
    void linearSolve(BoundaryType boundary, Grid& x, const Grid& x0, 
                    float a, float c, int iterations);
    
    // Utility methods
    void addSourceCircular(Grid& field, int i, int j, int radius, float value, float maxValue);
    bool isSolid(int i, int j) const;
};

} // namespace navier

#endif // FLUID_SOLVER_H
