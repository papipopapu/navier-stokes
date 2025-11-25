#ifndef GRID_H
#define GRID_H

#include <vector>
#include <stdexcept>
#include <cstring>
#include <algorithm>

namespace navier {

/**
 * @brief 2D grid container for fluid properties
 * Stores values in row-major order for cache efficiency
 */
class Grid {
private:
    size_t m_rows;    // Ni in original code
    size_t m_cols;    // Nj in original code
    std::vector<float> m_data;

public:
    Grid();
    Grid(size_t rows, size_t cols, float initialValue = 0.0f);
    
    // Accessors
    size_t rows() const { return m_rows; }
    size_t cols() const { return m_cols; }
    size_t size() const { return m_data.size(); }
    
    // Element access
    float& operator()(size_t i, size_t j);
    float operator()(size_t i, size_t j) const;
    
    // Safe element access with bounds checking
    float& at(size_t i, size_t j);
    float at(size_t i, size_t j) const;
    
    // Direct data access
    std::vector<float>& data() { return m_data; }
    const std::vector<float>& data() const { return m_data; }
    
    // Utility methods
    void resize(size_t rows, size_t cols, float initialValue = 0.0f);
    void fill(float value);
    void clear() { fill(0.0f); }
    
    // Drawing methods
    void drawLine(int i0, int j0, int length, bool isVertical, int thickness, float value);
    void drawRectangle(int i0, int j0, int width, int height, float value);
    void drawCircle(int centerI, int centerJ, int radius, float value);
    
    // Statistical methods
    float max() const;
    float min() const;
    float average() const;
    
    // Swap contents with another grid
    void swap(Grid& other);

private:
    inline size_t index(size_t i, size_t j) const {
        return i * m_cols + j;
    }
    
    inline bool isValidIndex(size_t i, size_t j) const {
        return i < m_rows && j < m_cols;
    }
};

/**
 * @brief Stores coordinates of solid cells for efficient iteration
 */
class SolidCoordinates {
private:
    std::vector<int> m_rows;
    std::vector<int> m_cols;

public:
    SolidCoordinates() = default;
    
    // Build from a solid map (non-zero cells are solid)
    void buildFromMap(const Grid& solidMap);
    
    // Accessors
    const std::vector<int>& rows() const { return m_rows; }
    const std::vector<int>& cols() const { return m_cols; }
    
    size_t count() const { return m_rows.size(); }
    bool empty() const { return m_rows.empty(); }
    
    void clear() {
        m_rows.clear();
        m_cols.clear();
    }
    
    // Get coordinate pair at index
    std::pair<int, int> at(size_t index) const {
        return {m_rows[index], m_cols[index]};
    }
};

/**
 * @brief Represents a line obstacle in the fluid domain
 */
struct Line {
    int i0, j0;           // Starting position
    int length;           // Length of the line
    int thickness;        // Thickness of the line
    bool isVertical;      // Direction of the line
    
    Line() : i0(0), j0(0), length(0), thickness(1), isVertical(false) {}
    
    Line(int i0, int j0, int length, bool isVertical, int thickness = 1)
        : i0(i0), j0(j0), length(length), thickness(thickness), isVertical(isVertical) {}
};

} // namespace navier

#endif // GRID_H
