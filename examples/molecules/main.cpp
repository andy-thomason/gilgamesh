
#include <meshutils/pdb_file.hpp>
#include <meshutils/mesh.hpp>
#include <meshutils/encoders/ply_encoder.hpp>
#include <glm/glm.hpp>
#include <ostream>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "2ptc.h"

#undef min

inline std::ostream & operator<<(std::ostream &os, const glm::vec3 &v) {
  return os << "vec3(" << v.x << ", " << v.y << ", " << v.z << ")";
}

class kdtree {
  struct sorter_t;
public:
  kdtree(const std::vector<glm::vec3> &pos) {
    std::vector<sorter_t> sorter;
    sorter.reserve(pos.size());
    
    for (auto &p : pos) {
      sorter.emplace_back(p, (int)(&p - pos.data()));
    }

    root_ = make_node(sorter.data(), 0, sorter.size(), 0);
  }

  int make_node(sorter_t *sorter, size_t begin, size_t end, int axis) {
    std::sort(
      sorter + begin, sorter + end,
      [axis, sorter](const sorter_t &a, const sorter_t &b) {
        return a.pos[axis] < b.pos[axis];
      }
    );

    size_t size = end - begin;
    int next_axis = axis+1 == 3 ? 0 : axis+1;
    size_t mid = begin + (size + 1) / 2;
    float split = mid == begin ? sorter[mid].pos[axis] : (sorter[mid-1].pos[axis] + sorter[mid].pos[axis]) * 0.5f;
    int left = 0;
    int right = 0;

    if (mid - begin >= 2) {
      left = make_node(sorter, begin, mid, next_axis);
    } else {
      assert(mid != begin);
      left = -1 - sorter[begin].index;
    }

    if (end - mid >= 2) {
      right = make_node(sorter, mid, end, next_axis);
    } else {
      assert(mid != end);
      right = -1 - sorter[mid].index;
    }

    int result = (int)nodes_.size();
    node_t new_node = { split, left, right };
    nodes_.push_back(new_node);

    //printf("%3d: %8.3f %3d %3d\n", result, split, left, right);
    return result;
  }

  // note that this is not thread safe
  void find_in_range(std::vector<size_t> &result, const glm::vec3 &min, const glm::vec3 &max) {
    result.resize(0);
    find(result, root_, min, max, 0);
  }

private:
  void find(std::vector<size_t> &result, int idx, const glm::vec3 &min, const glm::vec3 &max, int axis) {
    int next_axis = axis+1 == 3 ? 0 : axis+1;
    node_t &node = nodes_[idx];
    float split = node.split;
    //printf("f%3d: %8.3f %3d %3d %f..%f\n", idx, split, node.left, node.right, min[axis], max[axis]);
    if (min[axis] <= split) {
      if (node.left >= 0) {
        find(result, node.left, min, max, next_axis);
      } else {
        //printf("out: %d\n", -1-node.left);
        result.push_back(-1-node.left);
      }
    }
    if (max[axis] >= split) {
      if (node.right >= 0) {
        find(result, node.right, min, max, next_axis);
      } else {
        //printf("out: %d\n", -1-node.right);
        result.push_back(-1-node.right);
      }
    }
  }

  struct sorter_t {
    glm::vec3 pos;
    int index;

    sorter_t(glm::vec3 pos, int index) : pos(pos), index(index) {}
  };

  struct node_t {
    float split;
    int left;
    int right;
  };

  std::vector<node_t> nodes_;
  int root_;
};

int main() {
  meshutils::pdb_file pdb(
    (const uint8_t*)__2ptc_pdb, (const uint8_t*)__2ptc_pdb + sizeof(__2ptc_pdb)
  );

  std::string chains = pdb.chains();
  for (char chainID : chains) {
    std::vector<glm::vec3> pos = pdb.pos(chainID);
    std::vector<float> radii = pdb.radii(chainID);
    pos.resize(16);
    radii.resize(16);

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

    float water_radius = 1.4f;
    float grid_spacing = max_radius * 2;
    float recip_gs = 1.0f / grid_spacing;

    struct sphere {
      int x, y, z;
      int idx;
    };

    std::vector<sphere> spheres;
    for (size_t idx = 0; idx != pos.size(); ++idx) {
      glm::vec3 mins = (pos[idx] - min) - (radii[idx] + water_radius);
      glm::vec3 maxs = (pos[idx] - min) + (radii[idx] + water_radius);
      glm::ivec3 low = glm::ivec3(mins * recip_gs);
      glm::ivec3 high = glm::ivec3(maxs * recip_gs);
      for (int z = low.z; z <= high.z; ++z) {
        for (int y = low.y; y <= high.y; ++y) {
          for (int x = low.x; x <= high.x; ++x) {
            sphere s = { x, y, z, (int)idx };
            spheres.push_back(s);
          }
        }
      }
    }

    std::sort(
      spheres.begin(), spheres.end(),
      [](const sphere &a, const sphere &b) {
        return a.x < b.x + (a.y < b.y + (a.z < b.z));
      }
    );

    std::vector<std::pair<int, int>> pairs;
    int x = spheres[0].x, y = spheres[0].y, z = spheres[0].z;
    for (size_t i = 0; i != spheres.size(); ) {
      size_t j = i;
      for (; j != spheres.size() && spheres[i].x == spheres[j].x && spheres[i].y == spheres[j].y && spheres[i].z == spheres[j].z; ++j) {
        printf("%d,%d,%d,%d\n", spheres[j].x, spheres[j].y, spheres[j].z, spheres[j].idx);
      }
      pairs.emplace_back(spheres[i].idx, spheres[j].idx);
      printf("---\n");
      i = j;
    }

    for (size_t i = 0; i != pos.size(); ++i) {
      std::vector <size_t> indices;
      indices.resize(0);
      for (size_t j = i+1; j != pos.size(); ++j) {
        float min_distance = (radii[i] + water_radius*2) + radii[j];
        float d2 = glm::dot((pos[i] - pos[j]), (pos[i] - pos[j]));
        if (d2 * d2 <= min_distance * min_distance) {
          printf("%d, %d\n", i, j);
        }
      }

      //printf("[%3d] ", (int)i); for (auto x : indices) printf("%3d ", (int)x); printf("\n");
    }
  }

  meshutils::simple_mesh mesh;

  mesh.addCube();

  std::ofstream of("mesh.ply");

  meshutils::ply_encoder encoder;
  encoder.encode(mesh, of);
}
