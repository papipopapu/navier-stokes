#include "FluidSolver.h"
#include <algorithm>
#include <cmath>

namespace navier {

FluidSolver::FluidSolver(size_t rows, size_t cols)
    : m_rows(rows)
    , m_cols(cols)
    , m_time(0.0)
    , m_velocityU(rows, cols)
    , m_velocityV(rows, cols)
    , m_velocityU_prev(rows, cols)
    , m_velocityV_prev(rows, cols)
    , m_scalar(rows, cols)
    , m_scalar_prev(rows, cols)
    , m_pressure(rows, cols)
    , m_vorticity(rows, cols)
    , m_solidMap(rows, cols)
    , m_viscosity(SimulationConfig::DEFAULT_VISCOSITY)
    , m_diffusion(SimulationConfig::DEFAULT_DIFFUSION)
    , m_dissipation(SimulationConfig::DEFAULT_DISSIPATION)
    , m_boundaryU(BoundaryType::HORIZONTAL)
    , m_boundaryV(BoundaryType::VERTICAL)
    , m_boundaryScalar(BoundaryType::NONE)
    , m_boundaryPressure(BoundaryType::NONE)
{
}

void FluidSolver::setBoundaryConditions(BoundaryType u, BoundaryType v, 
                                        BoundaryType scalar, BoundaryType pressure) {
    m_boundaryU = u;
    m_boundaryV = v;
    m_boundaryScalar = scalar;
    m_boundaryPressure = pressure;
}

void FluidSolver::addObstacleLine(const Line& line) {
    m_solidMap.drawLine(line.i0, line.j0, line.length, line.isVertical, line.thickness, 1.0f);
    m_solidCoordinates.buildFromMap(m_solidMap);
}

void FluidSolver::addObstacleLines(std::initializer_list<Line> lines) {
    for (const auto& line : lines) {
        m_solidMap.drawLine(line.i0, line.j0, line.length, line.isVertical, line.thickness, 1.0f);
    }
    m_solidCoordinates.buildFromMap(m_solidMap);
}

void FluidSolver::clearObstacles() {
    m_solidMap.clear();
    m_solidCoordinates.clear();
}

void FluidSolver::addSourceCircular(Grid& field, int i, int j, int radius, float value, float maxValue) {
    int radiusSquared = radius * radius;
    for (int di = -radius; di <= radius; ++di) {
        for (int dj = -radius; dj <= radius; ++dj) {
            if (di * di + dj * dj <= radiusSquared) {
                int ni = i + di;
                int nj = j + dj;
                if (ni >= 0 && ni < static_cast<int>(m_rows) && 
                    nj >= 0 && nj < static_cast<int>(m_cols) && 
                    !isSolid(ni, nj)) {
                    field(ni, nj) = std::clamp(field(ni, nj) + value, -maxValue, maxValue);
                }
            }
        }
    }
}

void FluidSolver::addScalar(int i, int j, int radius, float value, float maxValue) {
    addSourceCircular(m_scalar, i, j, radius, value, maxValue);
}

void FluidSolver::addPressure(int i, int j, int radius, float value, float maxValue) {
    addSourceCircular(m_pressure, i, j, radius, value, maxValue);
}

void FluidSolver::addVelocity(int i, int j, int radius, float u, float v, float maxU, float maxV) {
    addSourceCircular(m_velocityU, i, j, radius, u, maxU);
    addSourceCircular(m_velocityV, i, j, radius, v, maxV);
}

bool FluidSolver::isSolid(int i, int j) const {
    if (i < 0 || i >= static_cast<int>(m_rows) || 
        j < 0 || j >= static_cast<int>(m_cols)) {
        return true;
    }
    return m_solidMap(i, j) != 0.0f;
}

void FluidSolver::step(float dt) {
    velocityStep(dt);
    scalarStep(dt);
    computeVorticity();
    m_time += dt;
}

void FluidSolver::velocityStep(float dt) {
    // Add viscous diffusion
    std::swap(m_velocityU, m_velocityU_prev);
    std::swap(m_velocityV, m_velocityV_prev);
    diffuse(m_boundaryU, dt, m_viscosity, m_velocityU, m_velocityU_prev);
    diffuse(m_boundaryV, dt, m_viscosity, m_velocityV, m_velocityV_prev);
    
    // Project to ensure incompressibility
    project();
    
    // Advect velocity by itself
    std::swap(m_velocityU, m_velocityU_prev);
    std::swap(m_velocityV, m_velocityV_prev);
    advect(m_boundaryU, dt, m_velocityU, m_velocityU_prev);
    advect(m_boundaryV, dt, m_velocityV, m_velocityV_prev);
    
    // Project again
    project();
    
    // Apply dissipation
    applyDissipation(m_velocityU, dt);
    applyDissipation(m_velocityV, dt);
}

void FluidSolver::scalarStep(float dt) {
    // Diffuse scalar
    std::swap(m_scalar, m_scalar_prev);
    diffuse(m_boundaryScalar, dt, m_diffusion, m_scalar, m_scalar_prev);
    
    // Advect scalar with velocity field
    std::swap(m_scalar, m_scalar_prev);
    advect(m_boundaryScalar, dt, m_scalar, m_scalar_prev);
    
    // Apply dissipation
    applyDissipation(m_scalar, dt);
}

void FluidSolver::advect(BoundaryType boundary, float dt, Grid& field, const Grid& field_prev) {
    float dt0 = dt * std::max(m_rows, m_cols);
    
    for (size_t i = 1; i < m_rows - 1; ++i) {
        for (size_t j = 1; j < m_cols - 1; ++j) {
            if (isSolid(i, j)) continue;
            
            // Backtrace along velocity field
            float x = static_cast<float>(i) - dt0 * m_velocityU_prev(i, j);
            float y = static_cast<float>(j) - dt0 * m_velocityV_prev(i, j);
            
            // Clamp to grid boundaries
            x = std::clamp(x, 0.5f, m_rows - 1.5f);
            y = std::clamp(y, 0.5f, m_cols - 1.5f);
            
            // Bilinear interpolation
            int i0 = static_cast<int>(x);
            int i1 = i0 + 1;
            int j0 = static_cast<int>(y);
            int j1 = j0 + 1;
            
            float s1 = x - i0;
            float s0 = 1.0f - s1;
            float t1 = y - j0;
            float t0 = 1.0f - t1;
            
            field(i, j) = s0 * (t0 * field_prev(i0, j0) + t1 * field_prev(i0, j1)) +
                         s1 * (t0 * field_prev(i1, j0) + t1 * field_prev(i1, j1));
        }
    }
    
    applyBoundaryConditions(boundary, field);
    applySolidBoundaries(boundary, field);
}

void FluidSolver::diffuse(BoundaryType boundary, float dt, float diffCoeff, 
                         Grid& field, const Grid& field_prev) {
    float a = dt * diffCoeff * m_rows * m_cols;
    linearSolve(boundary, field, field_prev, a, 1.0f + 4.0f * a, 
                SimulationConfig::GAUSS_SEIDEL_ITERATIONS);
}

void FluidSolver::linearSolve(BoundaryType boundary, Grid& x, const Grid& x0, 
                             float a, float c, int iterations) {
    float cRecip = 1.0f / c;
    
    for (int iter = 0; iter < iterations; ++iter) {
        for (size_t i = 1; i < m_rows - 1; ++i) {
            for (size_t j = 1; j < m_cols - 1; ++j) {
                if (isSolid(i, j)) continue;
                
                x(i, j) = (x0(i, j) + a * (x(i-1, j) + x(i+1, j) + 
                                           x(i, j-1) + x(i, j+1))) * cRecip;
            }
        }
        applyBoundaryConditions(boundary, x);
        applySolidBoundaries(boundary, x);
    }
}

void FluidSolver::project() {
    // Compute divergence
    Grid divergence(m_rows, m_cols);
    
    for (size_t i = 1; i < m_rows - 1; ++i) {
        for (size_t j = 1; j < m_cols - 1; ++j) {
            if (isSolid(i, j)) continue;
            
            divergence(i, j) = -0.5f * (m_velocityU(i+1, j) - m_velocityU(i-1, j) +
                                        m_velocityV(i, j+1) - m_velocityV(i, j-1)) / 
                                        std::max(m_rows, m_cols);
        }
    }
    applyBoundaryConditions(BoundaryType::NONE, divergence);
    applySolidBoundaries(BoundaryType::NONE, divergence);
    
    // Solve for pressure
    m_pressure.clear();
    linearSolve(m_boundaryPressure, m_pressure, divergence, 1.0f, 4.0f,
                SimulationConfig::GAUSS_SEIDEL_ITERATIONS);
    
    // Subtract pressure gradient from velocity
    for (size_t i = 1; i < m_rows - 1; ++i) {
        for (size_t j = 1; j < m_cols - 1; ++j) {
            if (isSolid(i, j)) continue;
            
            m_velocityU(i, j) -= 0.5f * (m_pressure(i+1, j) - m_pressure(i-1, j)) * 
                                std::max(m_rows, m_cols);
            m_velocityV(i, j) -= 0.5f * (m_pressure(i, j+1) - m_pressure(i, j-1)) * 
                                std::max(m_rows, m_cols);
        }
    }
    applyBoundaryConditions(m_boundaryU, m_velocityU);
    applyBoundaryConditions(m_boundaryV, m_velocityV);
    applySolidBoundaries(m_boundaryU, m_velocityU);
    applySolidBoundaries(m_boundaryV, m_velocityV);
}

void FluidSolver::applyDissipation(Grid& field, float dt) {
    float factor = std::exp(-m_dissipation * dt);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_cols; ++j) {
            field(i, j) *= factor;
        }
    }
}

void FluidSolver::applyBoundaryConditions(BoundaryType boundary, Grid& field) {
    // Set boundary cells based on type
    for (size_t i = 1; i < m_rows - 1; ++i) {
        field(i, 0) = (boundary == BoundaryType::VERTICAL) ? -field(i, 1) : field(i, 1);
        field(i, m_cols-1) = (boundary == BoundaryType::VERTICAL) ? -field(i, m_cols-2) : field(i, m_cols-2);
    }
    
    for (size_t j = 1; j < m_cols - 1; ++j) {
        field(0, j) = (boundary == BoundaryType::HORIZONTAL) ? -field(1, j) : field(1, j);
        field(m_rows-1, j) = (boundary == BoundaryType::HORIZONTAL) ? -field(m_rows-2, j) : field(m_rows-2, j);
    }
    
    // Corners
    field(0, 0) = 0.5f * (field(1, 0) + field(0, 1));
    field(0, m_cols-1) = 0.5f * (field(1, m_cols-1) + field(0, m_cols-2));
    field(m_rows-1, 0) = 0.5f * (field(m_rows-2, 0) + field(m_rows-1, 1));
    field(m_rows-1, m_cols-1) = 0.5f * (field(m_rows-2, m_cols-1) + field(m_rows-1, m_cols-2));
}

void FluidSolver::applySolidBoundaries(BoundaryType boundary, Grid& field) {
    // Apply boundary conditions at solid cells
    for (size_t idx = 0; idx < m_solidCoordinates.count(); ++idx) {
        auto [i, j] = m_solidCoordinates.at(idx);
        
        int count = 0;
        float sum = 0.0f;
        
        // Average from neighboring fluid cells
        if (i > 0 && !isSolid(i-1, j)) { sum += field(i-1, j); count++; }
        if (i < static_cast<int>(m_rows-1) && !isSolid(i+1, j)) { sum += field(i+1, j); count++; }
        if (j > 0 && !isSolid(i, j-1)) { sum += field(i, j-1); count++; }
        if (j < static_cast<int>(m_cols-1) && !isSolid(i, j+1)) { sum += field(i, j+1); count++; }
        
        if (count > 0) {
            field(i, j) = (boundary == BoundaryType::NONE) ? sum / count : -sum / count;
        }
    }
}

void FluidSolver::computeVorticity() {
    for (size_t i = 1; i < m_rows - 1; ++i) {
        for (size_t j = 1; j < m_cols - 1; ++j) {
            float dv_dx = (m_velocityV(i+1, j) - m_velocityV(i-1, j)) * 0.5f;
            float du_dy = (m_velocityU(i, j+1) - m_velocityU(i, j-1)) * 0.5f;
            m_vorticity(i, j) = dv_dx - du_dy;
        }
    }
}

void FluidSolver::applyPeriodicForcing(PeriodicFunction forceU, PeriodicFunction forceV,
                                      PeriodicFunction forceScalar, PeriodicFunction forcePressure, 
                                      float dt) {
    if (forceU) forceU(m_velocityU, m_time, dt);
    if (forceV) forceV(m_velocityV, m_time, dt);
    if (forceScalar) forceScalar(m_scalar, m_time, dt);
    if (forcePressure) forcePressure(m_pressure, m_time, dt);
}

} // namespace navier
