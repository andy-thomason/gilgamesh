////////////////////////////////////////////////////////////////////////////////
//
// (C) Andy Thomason 2016
//
// meshutils: box geometry class
// 

#ifndef MESHUTILS_BOX_INCLUDED
#define MESHUTILS_BOX_INCLUDED

#include <glm/glm.hpp>

namespace meshutils {

// An axis aligned box centred on the origin from -half_extent.x,y,z .. +half_extent.x,y,z
class box {
public:
  box(glm::vec3 half_extent = glm::vec3(1)) : half_extent_(half_extent) {
  }

  template <class Mesh>
  size_t build(Mesh &mesh, const glm::mat4 &transform = glm::mat4()) {
    size_t result = mesh.indices().size();

    //   6     7
    // 2 4  3  5
    // 0    1 
    static const int8_t cpos[8][3] = {
      { -1, -1, -1 },
      {  1, -1, -1 },
      { -1,  1, -1 },
      {  1,  1, -1 },
      { -1, -1,  1 },
      {  1, -1,  1 },
      { -1,  1,  1 },
      {  1,  1,  1 }
    };

    static const int8_t cnormal[6][3] = {
      { -1,  0,  0 },
      {  1,  0,  0 },
      {  0, -1,  0 },
      {  0,  1,  0 },
      {  0,  0, -1 },
      {  0,  0,  1 },
    };

    static const int8_t cuv[4][2] = {
      { 0, 0 },
      { 1, 0 },
      { 0, 1 },
      { 1, 1 },
    };

    //   6     7
    // 2 4  3  5
    // 0    1 
    static const uint8_t cube[] = {
      0, 4, 2, 6,
      1, 5, 3, 7,
      4, 5, 0, 1,
      2, 3, 6, 7,
      0, 1, 2, 3,
      5, 4, 7, 6,
    };

    for (int face = 0; face != 6; ++face) {
      glm::vec3 normal = glm::vec3(cnormal[face][0], cnormal[face][1], cnormal[face][2]);
      // four vertices per face
      for (int v = 0; v != 4; ++v) {
        size_t idx = cube[face*4+v];
        glm::vec3 pos = glm::vec3(cpos[idx][0], cpos[idx][1], cpos[idx][2]);
        glm::vec2 uv = glm::vec2(cuv[v][0], cuv[v][1]);
        mesh.addVertexTransformed(transform, pos * half_extent_, normal, uv);
      }

      // two triangles per face.
      mesh.addIndex(result + face*4 + 0);
      mesh.addIndex(result + face*4 + 2);
      mesh.addIndex(result + face*4 + 1);
      mesh.addIndex(result + face*4 + 1);
      mesh.addIndex(result + face*4 + 2);
      mesh.addIndex(result + face*4 + 3);
    }

    // result has 36 indices and 24 vertices
    return result;
  }
private:
  glm::vec3 half_extent_;
};

} // meshutils

#endif
