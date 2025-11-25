#include "Grid.h"
#include <cmath>
#include <limits>

namespace navier {

// Grid implementation
Grid::Grid() : m_rows(0), m_cols(0) {}

Grid::Grid(size_t rows, size_t cols, float initialValue)
    : m_rows(rows), m_cols(cols), m_data(rows * cols, initialValue) {}

float& Grid::operator()(size_t i, size_t j) {
    return m_data[index(i, j)];
}

float Grid::operator()(size_t i, size_t j) const {
    return m_data[index(i, j)];
}

float& Grid::at(size_t i, size_t j) {
    if (!isValidIndex(i, j)) {
        throw std::out_of_range("Grid index out of bounds");
    }
    return m_data[index(i, j)];
}

float Grid::at(size_t i, size_t j) const {
    if (!isValidIndex(i, j)) {
        throw std::out_of_range("Grid index out of bounds");
    }
    return m_data[index(i, j)];
}

void Grid::resize(size_t rows, size_t cols, float initialValue) {
    m_rows = rows;
    m_cols = cols;
    m_data.assign(rows * cols, initialValue);
}

void Grid::fill(float value) {
    std::fill(m_data.begin(), m_data.end(), value);
}

void Grid::drawLine(int i0, int j0, int length, bool isVertical, int thickness, float value) {
    if (isVertical) {
        for (int i = i0; i < i0 + length && i < static_cast<int>(m_rows); ++i) {
            for (int t = 0; t < thickness; ++t) {
                int j = j0 + t;
                if (j >= 0 && j < static_cast<int>(m_cols) && i >= 0) {
                    (*this)(i, j) = value;
                }
            }
        }
    } else {
        for (int j = j0; j < j0 + length && j < static_cast<int>(m_cols); ++j) {
            for (int t = 0; t < thickness; ++t) {
                int i = i0 + t;
                if (i >= 0 && i < static_cast<int>(m_rows) && j >= 0) {
                    (*this)(i, j) = value;
                }
            }
        }
    }
}

void Grid::drawRectangle(int i0, int j0, int width, int height, float value) {
    for (int i = i0; i < i0 + height && i < static_cast<int>(m_rows); ++i) {
        for (int j = j0; j < j0 + width && j < static_cast<int>(m_cols); ++j) {
            if (i >= 0 && j >= 0) {
                (*this)(i, j) = value;
            }
        }
    }
}

void Grid::drawCircle(int centerI, int centerJ, int radius, float value) {
    int radiusSquared = radius * radius;
    for (int i = centerI - radius; i <= centerI + radius; ++i) {
        for (int j = centerJ - radius; j <= centerJ + radius; ++j) {
            if (i >= 0 && i < static_cast<int>(m_rows) && 
                j >= 0 && j < static_cast<int>(m_cols)) {
                int di = i - centerI;
                int dj = j - centerJ;
                if (di * di + dj * dj <= radiusSquared) {
                    (*this)(i, j) = value;
                }
            }
        }
    }
}

float Grid::max() const {
    if (m_data.empty()) return 0.0f;
    return *std::max_element(m_data.begin(), m_data.end());
}

float Grid::min() const {
    if (m_data.empty()) return 0.0f;
    return *std::min_element(m_data.begin(), m_data.end());
}

float Grid::average() const {
    if (m_data.empty()) return 0.0f;
    float sum = 0.0f;
    for (float val : m_data) {
        sum += val;
    }
    return sum / m_data.size();
}

void Grid::swap(Grid& other) {
    std::swap(m_rows, other.m_rows);
    std::swap(m_cols, other.m_cols);
    m_data.swap(other.m_data);
}

// SolidCoordinates implementation
void SolidCoordinates::buildFromMap(const Grid& solidMap) {
    m_rows.clear();
    m_cols.clear();
    
    for (size_t i = 0; i < solidMap.rows(); ++i) {
        for (size_t j = 0; j < solidMap.cols(); ++j) {
            if (solidMap(i, j) != 0.0f) {
                m_rows.push_back(static_cast<int>(i));
                m_cols.push_back(static_cast<int>(j));
            }
        }
    }
}

} // namespace navier
