#ifndef RENDERER_H
#define RENDERER_H

#include "Grid.h"
#include "Config.h"
#include <SFML/Graphics.hpp>
#include <vector>
#include <memory>

namespace navier {

/**
 * @brief Color mapping strategies for visualization
 */
struct ColorMap {
    float r, g, b;
    uint8_t alpha;
    bool useGradient;
    bool normalize;
    
    ColorMap()
        : r(VisualizationConfig::DEFAULT_R)
        , g(VisualizationConfig::DEFAULT_G)
        , b(VisualizationConfig::DEFAULT_B)
        , alpha(VisualizationConfig::DEFAULT_ALPHA)
        , useGradient(true)
        , normalize(true)
    {}
};

/**
 * @brief Handles all rendering and visualization of fluid simulation
 * Separates rendering logic from simulation logic
 */
class Renderer {
private:
    size_t m_rows;
    size_t m_cols;
    float m_screenScale;
    
    std::unique_ptr<sf::RenderWindow> m_window;
    std::vector<sf::Uint8> m_pixels;
    sf::Texture m_texture;
    sf::Sprite m_sprite;
    
    ColorMap m_colorMap;
    VisualizationConfig::DrawMode m_drawMode;
    float m_drawMin;
    float m_drawMax;

public:
    Renderer(size_t rows, size_t cols, float screenScale = 1.0f);
    ~Renderer() = default;
    
    // Window management
    bool isOpen() const { return m_window && m_window->isOpen(); }
    void close() { if (m_window) m_window->close(); }
    sf::RenderWindow& window() { return *m_window; }
    
    // Color map configuration
    void setColorMap(const ColorMap& colorMap) { m_colorMap = colorMap; }
    void setColor(float r, float g, float b, uint8_t alpha = 255);
    void setGradientMode(bool useGradient) { m_colorMap.useGradient = useGradient; }
    void setNormalization(bool normalize) { m_colorMap.normalize = normalize; }
    
    // Drawing mode and range
    void setDrawMode(VisualizationConfig::DrawMode mode) { m_drawMode = mode; }
    void setDrawRange(float minVal, float maxVal);
    VisualizationConfig::DrawMode getDrawMode() const { return m_drawMode; }
    
    // Rendering methods
    void render(const Grid& field, const Grid& solidMap);
    void renderWithOptional(const Grid& field, const Grid& optional, const Grid& solidMap);
    void display();
    void clear();
    
    // Utility
    sf::Vector2i getMousePosition() const;
    void gridToScreen(int i, int j, int& x, int& y) const;
    void screenToGrid(int x, int y, int& i, int& j) const;

private:
    void updatePixels(const Grid& field, const Grid& solidMap, const Grid* optional = nullptr);
    sf::Color mapValueToColor(float value, float min, float max) const;
    sf::Color getGradientColor(float normalized) const;
    sf::Color getSolidColor(float value) const;
};

} // namespace navier

#endif // RENDERER_H
