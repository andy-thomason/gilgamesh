////////////////////////////////////////////////////////////////////////////////
//
// (C) Andy Thomason 2016
//
// Example using the basic_mesh class to genererate solvent excluded
// surfaces for molecules.
//


#include <meshutils/pdb_file.hpp>
#include <meshutils/mesh.hpp>
#include <meshutils/encoders/ply_encoder.hpp>

#include <glm/glm.hpp>

#include <ostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <thread>
#include <future>
#include <numeric>


#undef min

inline std::ostream & operator<<(std::ostream &os, const glm::vec3 &v) {
  return os << "vec3(" << v.x << ", " << v.y << ", " << v.z << ")";
}

template <class... P>
const char *fmt(const char *str, P... params) {
  static char buf[256];
  snprintf(buf, sizeof(buf), str, params...);
  return buf;
}

template <class F>
void par_for(int begin, int end, F fn) {
  std::atomic<int> idx;
  idx = begin;
  int num_cpus = std::thread::hardware_concurrency();
  std::vector<std::future<void>> futures(num_cpus);
  for (int cpu = 0; cpu != num_cpus; ++cpu) {
    futures[cpu] = std::async(std::launch::async, [cpu, &idx, end, &fn]() {
      for (;;) {
        int i = idx++;
        printf("[%d %d]", i, cpu);
        fflush(stdout);
        if (i >= end) break;
        fn(i);
      }
    });
  }
  for (int cpu = 0; cpu != num_cpus; ++cpu) {
    futures[cpu].get();
  }
  printf("\n");
}

int main() {
  std::ifstream file(CMAKE_SOURCE "/examples/data/2PTC.pdb", std::ios_base::binary);
  std::vector<uint8_t> text;
  if (!file.eof() && !file.fail()) {
    file.seekg(0, std::ios_base::end);
    text.resize((size_t)file.tellg());

    file.seekg(0, std::ios_base::beg);
    file.read((char*)text.data(), text.size());
  } 
  
  meshutils::pdb_file pdb(text.data(), text.data() + text.size());

  meshutils::ply_encoder encoder;
  
  std::string chains = pdb.chains();
  for (char chainID : chains) {
    std::vector<glm::vec3> pos = pdb.pos(chainID);
    std::vector<float> radii = pdb.radii(chainID);
    std::vector<glm::vec4> colors = pdb.colorsByFunction(chainID);

    // adjust centre of gravity
    glm::vec3 sum = std::accumulate(pos.begin(), pos.end(), glm::vec3(0, 0, 0));
    glm::vec3 cofg = sum * (1.0f/pos.size());
    for (auto &p : pos) {
      p -= cofg;
    }

    glm::vec3 min = pos[0];
    glm::vec3 max = pos[0];
    for (size_t i = 0; i != pos.size(); ++i) {
      glm::vec3 &p = pos[i];
      min = glm::min(min, p - glm::vec3(radii[i]));
      max = glm::max(max, p + glm::vec3(radii[i]));
    }

    float max_radius = 0;
    for (auto r : radii) {
      max_radius = std::max(max_radius, r);
    }

    float water_radius = 2.0f;
    float grid_spacing = 0.25;
    float recip_gs = 1.0f / grid_spacing;

    min -= water_radius + max_radius;
    max += water_radius + max_radius;

    int xdim = int((max.x - min.x) * recip_gs + 1);
    int ydim = int((max.y - min.y) * recip_gs + 1);
    int zdim = int((max.z - min.z) * recip_gs + 1);
    printf("%d x %d x %d\n", xdim, ydim, zdim);

    auto idx = [xdim, ydim](int x, int y, int z) {
      return ((z * (ydim+1)) + y) * (xdim+1) + x;
    };

    std::vector<float> accessible((xdim+1)*(ydim+1)*(zdim+1));
    par_for(0, zdim+1, [&](int z) {
      for (int y = 0; y != ydim+1; ++y) {
        for (int x = 0; x != xdim+1; ++x) {
          glm::vec3 xyz(x * grid_spacing + min.x, y * grid_spacing + min.y, z * grid_spacing + min.z);
          float value = 1e37f;
          for (size_t i = 0; i != pos.size(); ++i) {
            glm::vec3 p = pos[i];
            float r = radii[i];
            float d2 = glm::dot(xyz - p, xyz - p);
            float r2 = (r + water_radius) * (r + water_radius);
            value = std::min(value, d2 - r2);
          }
          accessible[idx(x, y, z)] = value;
        }
      }
    });

    auto fn = [&accessible, idx](int x, int y, int z) {
      return accessible[idx(x, y, z)];
    };

    auto gen = [&accessible, grid_spacing, min, idx](float x, float y, float z) {
      glm::vec3 xyz(x * grid_spacing + min.x, y * grid_spacing + min.y, z * grid_spacing + min.z);
      return meshutils::pos_mesh::vertex_t(xyz);
    };

    meshutils::pos_mesh amesh(xdim, ydim, zdim, fn, gen);
    std::ofstream of(fmt("accessible_%c.ply", chainID));
    encoder.encode(amesh, of);

    std::vector<glm::vec3> apos;
    auto *avertices = amesh.vertices();
    for (size_t i = 0; i != amesh.numVertices(); ++i) {
      apos.push_back(avertices[i].pos());
    }
    std::sort(
      apos.begin(), apos.end(),
      [](const glm::vec3 &a, const glm::vec3 &b) {
        return a.z < b.z;
      }
    );

    std::vector<float> excluded((xdim+1)*(ydim+1)*(zdim+1));
    float outside_value = -(water_radius * water_radius);
    par_for(0, zdim+1, [&](int z) {
      for (int y = 0; y != ydim+1; ++y) {
        for (int x = 0; x != xdim+1; ++x) {
          glm::vec3 xyz(x * grid_spacing + min.x, y * grid_spacing + min.y, z * grid_spacing + min.z);
          float value = 1e37f;
          if (accessible[idx(x, y, z)] < 0) {
            for (size_t i = 0; i != apos.size(); ++i) {
              glm::vec3 &p = apos[i];
              float d2 = glm::dot(xyz - p, xyz - p);
              float r2 = (water_radius) * (water_radius);
              value = std::min(value, d2 - r2);
            }
          } else {
            value = outside_value;
          }
          excluded[idx(x, y, z)] = value;
          if (z == 30 && y == 50) printf("%d %d %d %f\n", x, y, z, value);
        }
      }
    });

    auto efn = [&excluded, idx](int x, int y, int z) {
      return excluded[idx(x, y, z)];
    };

    auto egen = [&excluded, grid_spacing, min, idx](float x, float y, float z) {
      glm::vec3 xyz(x * grid_spacing + min.x, y * grid_spacing + min.y, z * grid_spacing + min.z);
      glm::vec3 normal(1, 0, 0);
      glm::vec2 uv(0, 0);
      glm::vec4 color = glm::vec4(1, 1, 1, 1);
      return meshutils::color_mesh::vertex_t(xyz, normal, uv, color);
    };

    meshutils::color_mesh emesh(xdim, ydim, zdim, efn, egen);

    std::ofstream eof(fmt("excluded_%c.ply", chainID));

    encoder.encode(emesh, eof);
  }
}
