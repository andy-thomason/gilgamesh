////////////////////////////////////////////////////////////////////////////////
//
// (C) Andy Thomason 2016
//
// meshutils: sphere geometry class
// 

#ifndef MESHUTILS_SPHERE_INCLUDED
#define MESHUTILS_SPHERE_INCLUDED

#include <glm/glm.hpp>

#include <random>

namespace meshutils {

// An axis aligned sphere centred on the origin from -half_extent.x,y,z .. +half_extent.x,y,z
class sphere {
public:
  sphere(float radius=1) : radius_(radius) {
  }

  template <class Mesh>
  size_t build(Mesh &mesh, const glm::mat4 &transform = glm::mat4(), int max_vertices=100) {
    size_t result = mesh.indices().size();

    int num_lattitude = 10;
    float length = 3.141592653589793f / (float)num_lattitude;

    for (int i = 0; i <= num_lattitude; ++i) {
      float angle = i * (3.141592653589793f / num_lattitude);
      float r = std::sin(angle);
      int num_longditude = std::max(1, int(r * (3.141592653589793f*2) / length));
      for (int j = 0; j != num_longditude; ++j) {
      }
    }


    /*for (int face = 0; face != 6; ++face) {
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
    }*/

    // result has 36 indices and 24 vertices
    return result;
  }
private:
  float radius_;
};

} // meshutils

#endif
