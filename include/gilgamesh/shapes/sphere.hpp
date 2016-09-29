////////////////////////////////////////////////////////////////////////////////
//
// (C) Andy Thomason 2016
//
// gilgamesh: sphere geometry class
// 

#ifndef MESHUTILS_SPHERE_INCLUDED
#define MESHUTILS_SPHERE_INCLUDED

#include <glm/glm.hpp>

namespace gilgamesh {

// A sphere centred on the origin.
class sphere {
public:
  sphere(float radius=1) : radius_(radius) {
  }

  // call this function to make a mesh
  // This works with the gilgamesh mesh, but is still generic.
  template <class Mesh>
  size_t build(Mesh &mesh, const glm::mat4 &transform = glm::mat4(), const glm::vec4 &color=glm::vec4(1), int num_lattitude=10, bool gen_uvs=false) {
    auto vertex = [&mesh, &transform, &color](const glm::vec3 &pos, const glm::vec3 &normal, const glm::vec2 &uv) {
      mesh.addVertexTransformed(transform, pos, normal, uv, color);
    };

    auto index = [&mesh](size_t idx) {
      mesh.addIndex(idx);
    };

    size_t first_index = mesh.vertices().size();

    buildMesh(vertex, index, first_index, num_lattitude, gen_uvs);
    
    // return the first index of the mesh
    return first_index;
  }

  // call this function to generate vertices and indices.
  template <class Vertex, class Index>
  void buildMesh(Vertex vertex, Index index, size_t first_index, int num_lattitude, bool gen_uvs) {
    // generate vertices
    float length = 3.141592653589793f / (float)num_lattitude;
    size_t pidx = first_index;
    int min_points = gen_uvs ? 6 : 1;
    for (int i = 0; i <= num_lattitude; ++i) {
      float phi = i * (3.141592653589793f / num_lattitude);
      float cosphi = std::cos(phi);
      float sinphi = std::sin(phi);
      int n = std::max(min_points, int(sinphi * (3.141592653589793f*2 / length)));
      int max_j = gen_uvs ? n + 1 : n;
      for (int j = 0; j < max_j; ++j) {
        float theta = j * (3.141592653589793f*2/n);
        float costheta = std::cos(theta);
        float sintheta = std::sin(theta);
        glm::vec3 pos(sinphi * costheta, cosphi, sinphi * sintheta);
        glm::vec2 uv(j * (1.0f/n), i * (1.0f / num_lattitude));
        vertex(pos * radius_, pos, uv);
      }
    }

    // generate indices
    float pphi = 0;
    int pn = 6;
    for (int i = 1; i <= num_lattitude; ++i) {
      float phi = i * (3.141592653589793f / num_lattitude);
      float sinphi = std::sin(phi);
      int n = std::max(min_points, int(sinphi * (3.141592653589793f*2 / length)));
      size_t idx = pidx + pn + (gen_uvs ? 1 : 0);

      // idx is the first vertex on the current row, pidx is the first on the previous row.
      if (gen_uvs) {
        for (int j = 0, pj = 0; j <= n && pj <= pn;) {
          index(idx + j);
          index(pidx + pj);
          if (j * pn < pj * n) { // equivalent to j / n < pj / pn
            index(idx + j + 1);
            ++j;
          } else {
            index(pidx + pj + 1);
            ++pj;
          }
        }
      } else {
        // Wrap at the meridian if not generating uv coordinates.
        // ie. there is one less point 
        auto wrap = [](int j, int n) {
          return j == n ? 0 : j;
        };

        for (int j = 0, pj = 0; j <= n && pj <= pn;) {
          index(idx + wrap(j, n));
          index(pidx + wrap(pj, pn));
          if (j * pn < pj * n) { // equivalent to j / n < pj / pn
            index(idx + wrap(j + 1, n));
            ++j;
          } else {
            index(pidx + wrap(pj + 1, pn));
            ++pj;
          }
        }
      }
      pidx = idx; pphi = phi; pn = n;
    }
  }
private:
  float radius_;
};

} // gilgamesh

#endif
