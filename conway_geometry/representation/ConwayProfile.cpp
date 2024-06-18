


#include "ConwayProfile.h"

#include <glm/glm.hpp>
#include <vector>

#include "../operations/geometryutils.h"

#include <sstream>
#include "SVGContext.h"

namespace conway::geometry {


std::string IfcProfile::DumpToOBJ() const {
  
  std::stringstream obj; 

  for ( const glm::dvec3& t : curve.points ) {
    
    obj << "v " << t.x << " " << t.y << " " << t.z << "\n";
  }

  size_t baseIndex = 1;

  for ( const IfcCurve& hole : holes ) {

    for ( const glm::dvec3& t : curve.points ) {
    
      obj << "v " << t.x << " " << t.y << " " << t.z << "\n";
    }
  }

  if ( curve.indices.size() > 1 ) {

      obj << "l";

      for ( uint16_t indice : curve.indices ) {

          obj << " " << indice + baseIndex;
      }

      obj << "\n";

      baseIndex += curve.points.size();
  
  } else if ( curve.points.size() > 1 ) {

      obj << "l";

      for ( size_t where = baseIndex, end = curve.points.size() + baseIndex; where < end; ++where ) {

          obj << " " << where;
      }
      
      obj << "\n";
  }

  for ( const IfcCurve& hole : holes ) {
    if ( hole.indices.size() > 1 ) {

      obj << "l";

      for ( uint16_t indice : curve.indices ) {

          obj << " " << indice + baseIndex;
      }

      obj << "\n";

      baseIndex += hole.points.size();
    
    } else if ( hole.points.size() > 1 ) {

      obj << "l";

      for ( size_t where = baseIndex, end = hole.points.size() + baseIndex; where < end; ++where ) {

          obj << " " << where;
      }
      
      obj << "\n";
    }
  }

  return obj.str();
}

std::string IfcProfile::DumpToSVG( const glm::dvec2& size, const glm::dvec2& offset ) const {
  
  glm::dvec2 min( FLT_MAX );
  glm::dvec2 max( -FLT_MAX );

  for ( const glm::dvec3& point : curve.points ) {

    min = glm::min( min, glm::dvec2( point ) );
    max = glm::max( max, glm::dvec2( point ) );
  }


  for ( const IfcCurve& hole : holes ) {
    
    for ( const glm::dvec3& point : hole.points ) {

      min = glm::min( min, glm::dvec2( point ) );
      max = glm::max( max, glm::dvec2( point ) );
    }
  }

  SVGContext svg( size, offset, min, max );

  svg.header();

  if ( curve.points.size() == 1 ) {

    svg.point( curve.points[ 0 ] );
  
  } else if ( !curve.indices.empty() ) {

    if ( curve.indices.size() == 1 ) {
      
      svg.point( curve.points[ curve.indices[ 0 ] ] );

    } else {

      for ( size_t where = 0, end = curve.indices.size() - 1; where < end; ++where ) {

        svg.line( curve.points[ curve.indices[ where ] ], curve.points[ curve.indices[ where + 1 ] ] );
      }
    }
  
  } else {

      for ( size_t where = 0, end = curve.points.size() - 1; where < end; ++where ) {

        svg.line( curve.points[ where ], curve.points[ where + 1 ] );
      }
  }

  for ( const IfcCurve& hole : holes ) {
    
    if ( hole.points.size() == 1 ) {

      svg.point( hole.points[ 0 ] );
    
    } else if ( !hole.indices.empty() ) {

      if ( hole.indices.size() == 1 ) {
        
        svg.point( hole.points[ hole.indices[ 0 ] ] );

      } else {

        for ( size_t where = 0, end = hole.indices.size() - 1; where < end; ++where ) {

          svg.line( hole.points[ hole.indices[ where ] ], hole.points[ hole.indices[ where + 1 ] ] );
        }
      }
    
    } else {

        for ( size_t where = 0, end = hole.points.size() - 1; where < end; ++where ) {

          svg.line( hole.points[ where ], hole.points[ where + 1 ] );
        }
    }

  }

  svg.trailer();

  return svg.str();
}


std::string IfcProfile::getType() const {
  return type;
}

IfcCurve IfcProfile::getCurve() const {
  return curve;
}

std::vector<IfcCurve> IfcProfile::getHoles() const {
  return holes;
}

bool IfcProfile::getIsConvex() const {
  return isConvex;
}

bool IfcProfile::getIsComposite() const {
  return isComposite;
}

std::vector<IfcProfile> IfcProfile::getProfiles() const {
  return profiles;
}

} // namespace conway::geometry