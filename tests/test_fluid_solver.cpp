#include "../core/FluidSolver.h"
#include <cassert>
#include <iostream>
#include <cmath>

using namespace navier;

void test_fluid_solver_construction() {
    std::cout << "Testing FluidSolver construction... ";
    
    FluidSolver solver(50, 60);
    
    assert(solver.rows() == 50);
    assert(solver.cols() == 60);
    assert(solver.time() == 0.0);
    
    // Check that grids are initialized
    const Grid& u = solver.velocityU();
    const Grid& v = solver.velocityV();
    const Grid& scalar = solver.scalar();
    const Grid& pressure = solver.pressure();
    const Grid& vorticity = solver.vorticity();
    const Grid& solidMap = solver.solidMap();
    
    assert(u.rows() == 50 && u.cols() == 60);
    assert(v.rows() == 50 && v.cols() == 60);
    assert(scalar.rows() == 50 && scalar.cols() == 60);
    assert(pressure.rows() == 50 && pressure.cols() == 60);
    assert(vorticity.rows() == 50 && vorticity.cols() == 60);
    assert(solidMap.rows() == 50 && solidMap.cols() == 60);
    
    std::cout << "PASSED\n";
}

void test_boundary_conditions() {
    std::cout << "Testing boundary conditions... ";
    
    FluidSolver solver(20, 20);
    
    // Test setting boundary conditions
    solver.setBoundaryConditions(
        BoundaryType::HORIZONTAL,
        BoundaryType::VERTICAL,
        BoundaryType::NONE,
        BoundaryType::NONE
    );
    
    // No crash means success - boundary conditions are applied during step
    
    std::cout << "PASSED\n";
}

void test_obstacle_management() {
    std::cout << "Testing obstacle management... ";
    
    FluidSolver solver(30, 30);
    
    // Add a single line obstacle
    Line line1(10, 10, 5, true, 1);
    solver.addObstacleLine(line1);
    
    // Check that obstacle is in solid map
    assert(solver.solidMap()(10, 10) != 0.0f);
    assert(solver.solidMap()(14, 10) != 0.0f);  // End of line
    assert(solver.solidMap()(15, 10) == 0.0f);  // Outside line
    
    // Clear obstacles
    solver.clearObstacles();
    assert(solver.solidMap()(10, 10) == 0.0f);
    
    // Add multiple lines at once
    solver.addObstacleLines({
        Line(5, 5, 10, true, 1),
        Line(5, 15, 10, true, 1)
    });
    
    assert(solver.solidMap()(5, 5) != 0.0f);
    assert(solver.solidMap()(5, 15) != 0.0f);
    
    std::cout << "PASSED\n";
}

void test_add_scalar() {
    std::cout << "Testing add scalar... ";
    
    FluidSolver solver(30, 30);
    
    // Add scalar at a point
    solver.addScalar(15, 15, 2, 100.0f, 1000.0f);
    
    // Check that scalar was added in the center
    assert(solver.scalar()(15, 15) > 0.0f);
    
    // Check that scalar respects max value
    for (int i = 0; i < 10; ++i) {
        solver.addScalar(15, 15, 2, 500.0f, 1000.0f);
    }
    assert(solver.scalar()(15, 15) <= 1000.0f);
    
    std::cout << "PASSED\n";
}

void test_add_velocity() {
    std::cout << "Testing add velocity... ";
    
    FluidSolver solver(30, 30);
    
    // Add velocity at a point
    solver.addVelocity(15, 15, 2, 50.0f, 30.0f, 100.0f, 100.0f);
    
    // Check that velocity was added
    assert(solver.velocityU()(15, 15) > 0.0f);
    assert(solver.velocityV()(15, 15) > 0.0f);
    
    std::cout << "PASSED\n";
}

void test_add_pressure() {
    std::cout << "Testing add pressure... ";
    
    FluidSolver solver(30, 30);
    
    // Add pressure at a point
    solver.addPressure(15, 15, 2, 1000.0f, 10000.0f);
    
    // Check that pressure was added
    assert(solver.pressure()(15, 15) > 0.0f);
    
    std::cout << "PASSED\n";
}

void test_simulation_step() {
    std::cout << "Testing simulation step... ";
    
    FluidSolver solver(30, 30);
    
    // Add some initial conditions
    solver.addScalar(15, 15, 3, 100.0f, 1000.0f);
    solver.addVelocity(15, 15, 3, 10.0f, 5.0f, 100.0f, 100.0f);
    
    // Initial time should be 0
    assert(solver.time() == 0.0);
    
    // Step the simulation
    float dt = 0.01f;
    solver.step(dt);
    
    // Time should have advanced
    assert(std::abs(solver.time() - dt) < 1e-6);
    
    // Multiple steps
    for (int i = 0; i < 10; ++i) {
        solver.step(dt);
    }
    
    // Time should have advanced by 11 * dt
    assert(std::abs(solver.time() - 11.0 * dt) < 1e-6);
    
    std::cout << "PASSED\n";
}

void test_vorticity_computation() {
    std::cout << "Testing vorticity computation... ";
    
    FluidSolver solver(30, 30);
    
    // Add some velocity to create vorticity
    // Add clockwise motion around center
    solver.addVelocity(14, 15, 1, 0.0f, -10.0f, 100.0f, 100.0f);
    solver.addVelocity(16, 15, 1, 0.0f, 10.0f, 100.0f, 100.0f);
    solver.addVelocity(15, 14, 1, 10.0f, 0.0f, 100.0f, 100.0f);
    solver.addVelocity(15, 16, 1, -10.0f, 0.0f, 100.0f, 100.0f);
    
    // Compute vorticity
    solver.computeVorticity();
    
    // There should be non-zero vorticity somewhere
    float maxVort = solver.vorticity().max();
    float minVort = solver.vorticity().min();
    assert(maxVort != 0.0f || minVort != 0.0f);
    
    std::cout << "PASSED\n";
}

void test_dissipation() {
    std::cout << "Testing dissipation... ";
    
    FluidSolver solver(30, 30);
    solver.setDissipation(0.5f);  // High dissipation for testing
    
    // Add scalar
    solver.addScalar(15, 15, 3, 100.0f, 1000.0f);
    float initialMax = solver.scalar().max();
    
    // Step multiple times
    for (int i = 0; i < 10; ++i) {
        solver.step(0.01f);
    }
    
    // Scalar should have decreased due to dissipation
    float finalMax = solver.scalar().max();
    assert(finalMax < initialMax);
    
    std::cout << "PASSED\n";
}

void test_obstacle_interaction() {
    std::cout << "Testing obstacle interaction... ";
    
    FluidSolver solver(30, 30);
    
    // Add an obstacle in the middle
    solver.addObstacleLine(Line(10, 15, 10, true, 2));
    
    // Add scalar near the obstacle
    solver.addScalar(10, 10, 2, 100.0f, 1000.0f);
    
    // Step the simulation
    for (int i = 0; i < 5; ++i) {
        solver.step(0.01f);
    }
    
    // Simulation should run without crashing
    // Obstacle cells should remain in solid map
    assert(solver.solidMap()(15, 15) != 0.0f);
    
    std::cout << "PASSED\n";
}

void test_physical_parameters() {
    std::cout << "Testing physical parameter setters... ";
    
    FluidSolver solver(20, 20);
    
    // Set various parameters (no getters to verify, but shouldn't crash)
    solver.setViscosity(0.5f);
    solver.setDiffusion(0.1f);
    solver.setDissipation(0.01f);
    
    // Step should work with different parameters
    solver.addScalar(10, 10, 2, 50.0f, 1000.0f);
    solver.step(0.01f);
    
    std::cout << "PASSED\n";
}

void test_periodic_forcing() {
    std::cout << "Testing periodic forcing... ";
    
    FluidSolver solver(20, 20);
    
    // Define a simple forcing function
    auto forceU = [](Grid& u, double /*time*/, double /*dt*/) {
        u(10, 10) += 1.0f;
    };
    
    // Apply periodic forcing
    solver.applyPeriodicForcing(forceU, nullptr, nullptr, nullptr, 0.01f);
    
    // Check that velocity was modified
    assert(solver.velocityU()(10, 10) == 1.0f);
    
    // Apply again
    solver.applyPeriodicForcing(forceU, nullptr, nullptr, nullptr, 0.01f);
    assert(solver.velocityU()(10, 10) == 2.0f);
    
    std::cout << "PASSED\n";
}

int main() {
    std::cout << "=== Running FluidSolver Tests ===\n";
    
    test_fluid_solver_construction();
    test_boundary_conditions();
    test_obstacle_management();
    test_add_scalar();
    test_add_velocity();
    test_add_pressure();
    test_simulation_step();
    test_vorticity_computation();
    test_dissipation();
    test_obstacle_interaction();
    test_physical_parameters();
    test_periodic_forcing();
    
    std::cout << "=== All FluidSolver tests passed! ===\n";
    return 0;
}
