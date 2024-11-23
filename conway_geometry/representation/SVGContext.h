#pragma once

#include <sstream>
#include <glm/glm.hpp>


namespace conway::geometry {

  struct SVGContext {

    std::stringstream stream;
    glm::dvec2 scale;
    glm::dvec2 bias;
    glm::dvec2 size;

    SVGContext(
      const glm::dvec2& _size,
      const glm::dvec2& offset,
      const glm::dvec2& min,
      const glm::dvec2& max ) :
        scale( ( _size - ( offset * 2.0 ) ) / ( max - min ) ),
        bias( ( offset / scale ) - min ),
        size( _size )
        {}

    template< typename input >
    SVGContext& operator<<( const input& value ) {
      stream << value;
      return *this;
    }

    SVGContext& header() {
      stream << "<svg width=\"" << size.x << "\" height=\"" << size.y << " \" xmlns=\"http://www.w3.org/2000/svg\">";
      return *this;
    }

    SVGContext& trailer() {
      stream << "</svg>";
      return *this;
    }

    glm::dvec2 toSvg( const glm::dvec2& from ) const {

      return ( from + bias ) * scale;
    }


    SVGContext& line( const glm::dvec3& p0, const glm::dvec3& p1, const char* color = "red", uint32_t width = 1 ) {
    
      return line( glm::dvec2( p0 ), glm::dvec2( p1 ), color, width );
    }
      

    SVGContext& line( const glm::dvec2& p0, const glm::dvec2& p1, const char* color = "red", uint32_t width = 1 ) {
      
      glm::dvec2 a = toSvg( p0 );
      glm::dvec2 b = toSvg( p1 );

      stream << "<line x1=\"" << a.x << "\" y1=\"" << a.y << "\" "
            << "x2=\"" << b.x << "\" y2=\"" << b.y << "\" "
            << "style = \"stroke:" << color << ";stroke-width:" << width << "\"/>";

      return *this;
    }

    SVGContext& point( const char* label, const glm::dvec2& p, const char* color = "blue", uint32_t width = 2 ) {
      
      glm::dvec2 a = toSvg( p );

      stream << "<circle cx = \"" << a.x << "\" cy = \""
            << a.y << "\" r = \"3\" style = \"stroke:"
            << color << ";stroke-width:" << width << "\" />"
            << "<text x=\"" << a.x + 2 * width << "\" y=\"" << a.y + width << "\">" << label << "</text>";

      return *this;
    }

    SVGContext& point( const glm::dvec2& p, const char* color = "blue", uint32_t width = 2 ) {
      
      glm::dvec2 a = toSvg( p );

      stream << "<circle cx = \"" << a.x << "\" cy = \""
            << a.y << "\" r = \"3\" style = \"stroke:"
            << color << ";stroke-width:" << width << "\" />";

      return *this;
    }

    std::string str() const {
      return stream.str();
    }

  };
}