#pragma once 

#include "structures/union_find.h"
#include <unordered_map>

#include <CDT/CDT.h>

namespace conway::geometry {

  class FaceCharter {
  public:  

    UnionFind< uint32_t > vertices;
    UnionFind< uint32_t > face_planes;

    std::unordered_map< std::pair< uint32_t, uint32_t >, uint32_t > edge_edge_vertices; 
    std::unordered_map< std::pair< uint32_t, uint32_t >, uint32_t > face_edge_vertices;

  private:


  };
}