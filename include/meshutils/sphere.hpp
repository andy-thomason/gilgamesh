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

    std::minstd_rand0 random;
    std::uniform_real_distribution<float> dis(-1, 1);

    std::vector<glm::vec3> vertices(max_vertices);
    for (int i = 0; i != max_vertices; ++i) {
      float x, y, z;
      do {
        x = dis(random);
        y = dis(random);
        z = dis(random);
      } while (x*x + y*y + z*z < 1e-6f);
      vertices[i] = glm::normalize(v);
    }

    for (int i = 0; i != max_vertices; ++i) {
      auto vi = vertices[i];
      float dmin = 1e37;
      int jmin = -1;
      for (int j = i + 1; j < max_vertices; ++j) {
        auto vj = vertices[j];
        auto dot = glm::dot(v0, vj);
        if (dot < dmin) {
          dmin = dot;
          jmin = j;
        }
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
