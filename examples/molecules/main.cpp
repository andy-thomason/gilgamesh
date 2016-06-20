
#include <meshutils/pdb_file.hpp>
#include <glm/glm.hpp>
#include <ostream>
#include <cmath>
#include <algorithm>

#include "2ptc.h"

#undef min

inline std::ostream & operator<<(std::ostream &os, const glm::vec3 &v) {
  return os << "vec3(" << v.x << ", " << v.y << ", " << v.z << ")";
}

class kdtree {
public:
  kdtree(const std::vector<glm::vec3> &pos) {
    std::vector<sorter_t> sorter;

    for (auto &p : pos) {
      sorter.emplace_back(p, &p - pos.data());
    }

    std::sort(
      sorter.data(), sorter.data() + sorter.size(),
      [](const sorter_t &a, const sorter_t &b) {
        return a.pos.x < b.pos.x;
      }
    );
  }

private:
  struct sorter_t {
    glm::vec3 pos;
    size_t index;

    sorter_t(glm::vec3 pos, size_t index) : pos(pos), index(index) {}
  };

  //std::vector<float> partitions;
};

int main() {
  meshutils::pdb_file pdb((const uint8_t*)__2ptc_pdb, (const uint8_t*)__2ptc_pdb + sizeof(__2ptc_pdb));

  std::string chains = pdb.chains();
  for (char chainID : chains) {
    std::vector<glm::vec3> pos = pdb.pos(chainID);
    std::vector<float> radii = pdb.radii(chainID);

    kdtree kd(pos);

    glm::vec3 min = pos[0];
    glm::vec3 max = pos[0];
    for (auto &p : pos) {
      min = glm::min(min, p);
      max = glm::max(max, p);
    }

    float spacing = 1.0f;
    std::cout << min << ", " << max << "\n";
    glm::ivec3 dim = glm::ivec3((max - min) / spacing) + glm::ivec3(1, 1, 1);

    int xdim = dim.x, ydim = dim.y, zdim = dim.z;
    std::vector<float> distance((xdim+1)*(ydim+1)*(zdim+1));

    float *value = distance.data();
    for (int z = 0; z <= zdim; ++z) {
      for (int y = 0; y <= ydim; ++y) {
        for (int x = 0; x <= xdim; ++x) {
          glm::vec3 p(x * spacing + min.x, y * spacing + min.y, z * spacing + min.z);
          float min_d2 = 1e30f;
          /*for (size_t i = 0; i != pos.size(); ++i) {
            float r = radii[i];
            float d2 = glm::dot(pos[i] - p, pos[i] - p) - r * r;
            min_d2 = std::min(min_d2, d2);
          }*/
          *value++ = min_d2;
        }
      }
    }
  }
}
