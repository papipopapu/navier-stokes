#include "../core/Grid.h"
#include <cassert>
#include <iostream>
#include <cmath>

using namespace navier;

void test_grid_construction() {
    std::cout << "Testing Grid construction... ";
    
    // Default constructor
    Grid g1;
    assert(g1.rows() == 0);
    assert(g1.cols() == 0);
    assert(g1.size() == 0);
    
    // Parameterized constructor
    Grid g2(10, 20);
    assert(g2.rows() == 10);
    assert(g2.cols() == 20);
    assert(g2.size() == 200);
    
    // Constructor with initial value
    Grid g3(5, 5, 1.5f);
    assert(g3.rows() == 5);
    assert(g3.cols() == 5);
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            assert(g3(i, j) == 1.5f);
        }
    }
    
    std::cout << "PASSED\n";
}

void test_grid_access() {
    std::cout << "Testing Grid element access... ";
    
    Grid g(10, 10);
    
    // Set and get values
    g(3, 5) = 42.0f;
    assert(g(3, 5) == 42.0f);
    
    g(0, 0) = 1.0f;
    assert(g(0, 0) == 1.0f);
    
    g(9, 9) = 99.0f;
    assert(g(9, 9) == 99.0f);
    
    // Test at() method
    g.at(5, 5) = 55.0f;
    assert(g.at(5, 5) == 55.0f);
    
    // Test out of bounds access
    try {
        g.at(10, 5);
        assert(false && "Should have thrown exception");
    } catch (const std::out_of_range&) {
        // Expected
    }
    
    std::cout << "PASSED\n";
}

void test_grid_utility_methods() {
    std::cout << "Testing Grid utility methods... ";
    
    Grid g(5, 5);
    
    // Test fill
    g.fill(3.0f);
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            assert(g(i, j) == 3.0f);
        }
    }
    
    // Test clear
    g.clear();
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            assert(g(i, j) == 0.0f);
        }
    }
    
    // Test resize
    g.resize(3, 4, 2.0f);
    assert(g.rows() == 3);
    assert(g.cols() == 4);
    assert(g.size() == 12);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            assert(g(i, j) == 2.0f);
        }
    }
    
    std::cout << "PASSED\n";
}

void test_grid_statistics() {
    std::cout << "Testing Grid statistics... ";
    
    Grid g(3, 3);
    g(0, 0) = 1.0f;
    g(0, 1) = 2.0f;
    g(0, 2) = 3.0f;
    g(1, 0) = 4.0f;
    g(1, 1) = 5.0f;
    g(1, 2) = 6.0f;
    g(2, 0) = 7.0f;
    g(2, 1) = 8.0f;
    g(2, 2) = 9.0f;
    
    assert(g.min() == 1.0f);
    assert(g.max() == 9.0f);
    assert(std::abs(g.average() - 5.0f) < 0.001f);
    
    std::cout << "PASSED\n";
}

void test_grid_swap() {
    std::cout << "Testing Grid swap... ";
    
    Grid g1(3, 3, 1.0f);
    Grid g2(5, 5, 2.0f);
    
    g1.swap(g2);
    
    assert(g1.rows() == 5);
    assert(g1.cols() == 5);
    assert(g1(0, 0) == 2.0f);
    
    assert(g2.rows() == 3);
    assert(g2.cols() == 3);
    assert(g2(0, 0) == 1.0f);
    
    std::cout << "PASSED\n";
}

void test_grid_drawing() {
    std::cout << "Testing Grid drawing methods... ";
    
    // Test drawLine
    Grid g(10, 10);
    g.drawLine(2, 2, 5, true, 1, 1.0f);  // Vertical line at column 2, from row 2, length 5
    assert(g(2, 2) == 1.0f);
    assert(g(3, 2) == 1.0f);
    assert(g(6, 2) == 1.0f);
    assert(g(7, 2) == 0.0f);  // Outside line
    
    g.clear();
    g.drawLine(2, 2, 5, false, 1, 1.0f);  // Horizontal line at row 2, from column 2, length 5
    assert(g(2, 2) == 1.0f);
    assert(g(2, 6) == 1.0f);
    assert(g(2, 7) == 0.0f);  // Outside line
    
    // Test drawRectangle
    g.clear();
    g.drawRectangle(3, 3, 2, 2, 5.0f);
    assert(g(3, 3) == 5.0f);
    assert(g(3, 4) == 5.0f);
    assert(g(4, 3) == 5.0f);
    assert(g(4, 4) == 5.0f);
    assert(g(5, 5) == 0.0f);  // Outside rectangle
    
    // Test drawCircle
    g.clear();
    g.drawCircle(5, 5, 2, 3.0f);
    assert(g(5, 5) == 3.0f);  // Center
    assert(g(5, 3) == 3.0f);  // Left edge
    assert(g(5, 7) == 3.0f);  // Right edge
    
    std::cout << "PASSED\n";
}

void test_solid_coordinates() {
    std::cout << "Testing SolidCoordinates... ";
    
    Grid solidMap(5, 5);
    solidMap(1, 1) = 1.0f;
    solidMap(2, 3) = 1.0f;
    solidMap(4, 4) = 1.0f;
    
    SolidCoordinates coords;
    coords.buildFromMap(solidMap);
    
    assert(coords.count() == 3);
    assert(!coords.empty());
    
    // Check the coordinates are correct
    auto [i0, j0] = coords.at(0);
    auto [i1, j1] = coords.at(1);
    auto [i2, j2] = coords.at(2);
    
    assert(i0 == 1 && j0 == 1);
    assert(i1 == 2 && j1 == 3);
    assert(i2 == 4 && j2 == 4);
    
    // Test clear
    coords.clear();
    assert(coords.count() == 0);
    assert(coords.empty());
    
    std::cout << "PASSED\n";
}

void test_line_struct() {
    std::cout << "Testing Line struct... ";
    
    // Default constructor
    Line l1;
    assert(l1.i0 == 0);
    assert(l1.j0 == 0);
    assert(l1.length == 0);
    assert(l1.thickness == 1);
    assert(l1.isVertical == false);
    
    // Parameterized constructor
    Line l2(10, 20, 30, true, 5);
    assert(l2.i0 == 10);
    assert(l2.j0 == 20);
    assert(l2.length == 30);
    assert(l2.isVertical == true);
    assert(l2.thickness == 5);
    
    std::cout << "PASSED\n";
}

int main() {
    std::cout << "=== Running Grid Tests ===\n";
    
    test_grid_construction();
    test_grid_access();
    test_grid_utility_methods();
    test_grid_statistics();
    test_grid_swap();
    test_grid_drawing();
    test_solid_coordinates();
    test_line_struct();
    
    std::cout << "=== All Grid tests passed! ===\n";
    return 0;
}
