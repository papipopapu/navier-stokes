#include "Renderer.h"
#include <algorithm>
#include <cmath>

namespace navier {

Renderer::Renderer(size_t rows, size_t cols, float screenScale)
    : m_rows(rows)
    , m_cols(cols)
    , m_screenScale(screenScale)
    , m_drawMode(VisualizationConfig::DrawMode::SCALAR)
    , m_drawMin(0.0f)
    , m_drawMax(VisualizationConfig::MAX_SCALAR)
{
    // Create window
    unsigned int windowWidth = static_cast<unsigned int>(cols * screenScale);
    unsigned int windowHeight = static_cast<unsigned int>(rows * screenScale);
    m_window = std::make_unique<sf::RenderWindow>(
        sf::VideoMode(windowWidth, windowHeight),
        "Navier-Stokes Fluid Simulation"
    );
    
    // Initialize pixel buffer
    m_pixels.resize(rows * cols * 4); // RGBA
    
    // Create texture and sprite
    m_texture.create(static_cast<unsigned int>(cols), static_cast<unsigned int>(rows));
    m_sprite.setTexture(m_texture);
    m_sprite.setScale(screenScale, screenScale);
}

void Renderer::setColor(float r, float g, float b, uint8_t alpha) {
    m_colorMap.r = r;
    m_colorMap.g = g;
    m_colorMap.b = b;
    m_colorMap.alpha = alpha;
}

void Renderer::setDrawRange(float minVal, float maxVal) {
    m_drawMin = minVal;
    m_drawMax = maxVal;
}

void Renderer::render(const Grid& field, const Grid& solidMap) {
    updatePixels(field, solidMap, nullptr);
    m_texture.update(m_pixels.data());
}

void Renderer::renderWithOptional(const Grid& field, const Grid& optional, const Grid& solidMap) {
    updatePixels(field, solidMap, &optional);
    m_texture.update(m_pixels.data());
}

void Renderer::display() {
    m_window->clear();
    m_window->draw(m_sprite);
    m_window->display();
}

void Renderer::clear() {
    m_window->clear();
}

sf::Vector2i Renderer::getMousePosition() const {
    return sf::Mouse::getPosition(*m_window);
}

void Renderer::gridToScreen(int i, int j, int& x, int& y) const {
    x = static_cast<int>(j * m_screenScale);
    y = static_cast<int>(i * m_screenScale);
}

void Renderer::screenToGrid(int x, int y, int& i, int& j) const {
    j = static_cast<int>(x / m_screenScale);
    i = static_cast<int>(y / m_screenScale);
}

void Renderer::updatePixels(const Grid& field, const Grid& solidMap, const Grid* optional) {
    float minVal = m_drawMin;
    float maxVal = m_drawMax;
    
    // Compute actual min/max if normalizing
    if (m_colorMap.normalize) {
        minVal = field.min();
        maxVal = field.max();
        if (std::abs(maxVal - minVal) < 1e-6f) {
            maxVal = minVal + 1.0f;
        }
    }
    
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_cols; ++j) {
            size_t pixelIndex = (i * m_cols + j) * 4;
            
            // Check if solid
            if (solidMap(i, j) != 0.0f) {
                // Draw solid as black
                m_pixels[pixelIndex] = 0;
                m_pixels[pixelIndex + 1] = 0;
                m_pixels[pixelIndex + 2] = 0;
                m_pixels[pixelIndex + 3] = 255;
            } else {
                float value = field(i, j);
                
                // Apply optional field overlay
                if (optional && (*optional)(i, j) != 0.0f) {
                    value += (*optional)(i, j);
                }
                
                sf::Color color = mapValueToColor(value, minVal, maxVal);
                m_pixels[pixelIndex] = color.r;
                m_pixels[pixelIndex + 1] = color.g;
                m_pixels[pixelIndex + 2] = color.b;
                m_pixels[pixelIndex + 3] = color.a;
            }
        }
    }
}

sf::Color Renderer::mapValueToColor(float value, float min, float max) const {
    // Ensure we have a valid range to avoid division by zero
    float range = max - min;
    if (std::abs(range) < 1e-6f) {
        max = min + 1.0f;
        range = 1.0f;
    }
    
    // Normalize value to [0, 1]
    float normalized = (value - min) / range;
    normalized = std::clamp(normalized, 0.0f, 1.0f);
    
    if (m_colorMap.useGradient) {
        return getGradientColor(normalized);
    } else {
        return getSolidColor(normalized);
    }
}

sf::Color Renderer::getGradientColor(float normalized) const {
    // Blue to white to red gradient
    sf::Color color;
    
    if (normalized < 0.5f) {
        // Blue to white
        float t = normalized * 2.0f;
        color.r = static_cast<uint8_t>(t * 255);
        color.g = static_cast<uint8_t>(t * 255);
        color.b = 255;
    } else {
        // White to red
        float t = (normalized - 0.5f) * 2.0f;
        color.r = 255;
        color.g = static_cast<uint8_t>((1.0f - t) * 255);
        color.b = static_cast<uint8_t>((1.0f - t) * 255);
    }
    
    color.a = m_colorMap.alpha;
    return color;
}

sf::Color Renderer::getSolidColor(float normalized) const {
    // Use the configured color with intensity based on normalized value
    return sf::Color(
        static_cast<uint8_t>(m_colorMap.r * normalized * 255),
        static_cast<uint8_t>(m_colorMap.g * normalized * 255),
        static_cast<uint8_t>(m_colorMap.b * normalized * 255),
        m_colorMap.alpha
    );
}

} // namespace navier
