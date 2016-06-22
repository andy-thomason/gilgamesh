
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
    sorter.reserve(pos.size());
    
    for (auto &p : pos) {
      sorter.emplace_back(p, &p - pos.data());
    }

    size_t size = sorter.size();
    auto *begin = sorter.data();
    auto *end = begin + size;

    int axis = 0;
    size_t step = size;
    for (;;) {
      for (size_t i = 0; i < size; i += step) {
        size_t max_j = std::min(i + step, size);
        for (size_t j = i; j < max_j; ++j) {
          begin[j].partition = i;
        }
      }

      std::sort(
        begin, end,
        [axis, begin](const sorter_t &a, const sorter_t &b) {
          return a.partition < (b.partition + (a.pos[axis] < b.pos[axis]));
        }
      );

      if (step == 1) {
        printf("---\n");
        for (auto *p = begin; p != end; ++p) {
          printf("%5d %5d %f %f %f %d\n", (int)(p - begin), (int)(p->partition), p->pos.x, p->pos.y, p->pos.z, (int)p->index);
        }
      }


      if (step == 1) break;
      
      for (size_t lower = 0; lower < size; lower += step) {
        size_t upper = std::min(lower + step, size);
        size_t mid = lower + (step + 1) / 2;
        partition_t p = { lower, upper, mid, begin[mid].pos[axis] };
        printf("p%d %d %d %d %f\n", (int)partitions_.size(), (int)lower, (int)mid, (int)upper, begin[mid].pos[axis]);
        partitions_.push_back(p);
      }
      axis = axis == 2 ? 0 : axis + 1;
      step = (step + 1) / 2;
    }

    in_starts_.reserve(size);
    out_starts_.reserve(size);
    size_ = size;
  }

  // note that this is not thread safe
  const std::vector <size_t> &find_in_range(const glm::vec3 &min, const glm::vec3 &max) {
    int axis = 0;
    size_t step = size_;
    in_starts_.resize(0);
    in_starts_.emplace_back(0);
    for (;;) {
      out_starts_.resize(0);
      float minx = min[axis];
      float maxx = max[axis];
      printf("s=%d %f..%f\n", (int)step, minx, maxx);
      for (size_t idx : in_starts_) {
        partition_t p = partitions_[idx];
        size_t lower = p.lower;
        size_t upper = p.upper;
        size_t mid = p.mid;
        float mid_value = p.mid_value;
        printf("%d..%d..%d mid_value=%f\n", (int)lower, (int)mid, (int)upper, mid_value);
        if (minx <= mid_value) {
          out_starts_.emplace_back(lower);
        }
        if (maxx >= mid_value) {
          out_starts_.emplace_back(mid);
        }
      }

      if (step == 1) break;

      axis = axis == 2 ? 0 : axis + 1;
      step = (step + 1) / 2;
      std::swap(in_starts_, out_starts_);
    }
    return out_starts_;
  }

private:
  struct sorter_t {
    glm::vec3 pos;
    size_t index;
    size_t partition;

    sorter_t(glm::vec3 pos, size_t index) : pos(pos), index(index) {}
  };

  struct partition_t {
    size_t lower;
    size_t upper;
    size_t mid;
    float mid_value;
  };

  std::vector<size_t> in_starts_;
  std::vector<size_t> out_starts_;
  std::vector<partition_t> partitions_;
  size_t size_;
};

int main() {
  meshutils::pdb_file pdb(
    (const uint8_t*)__2ptc_pdb, (const uint8_t*)__2ptc_pdb + sizeof(__2ptc_pdb)
  );

  std::string chains = pdb.chains();
  for (char chainID : chains) {
    std::vector<glm::vec3> pos = pdb.pos(chainID);
    std::vector<float> radii = pdb.radii(chainID);

    kdtree kd(pos);

    const std::vector <size_t> &res = kd.find_in_range(pos[0], pos[0]);
    return 0;

    glm::vec3 min = pos[0];
    glm::vec3 max = pos[0];
    for (auto &p : pos) {
      min = glm::min(min, p);
      max = glm::max(max, p);
    }

    float water_radius = 1.52f;
    size_t count = 0;
    size_t count2 = 0;
    for (size_t i = 0; i != pos.size(); ++i) {
      glm::vec3 min = pos[i] - glm::vec3(radii[i] + water_radius*2 + 3);
      glm::vec3 max = pos[i] + glm::vec3(radii[i] + water_radius*2 + 3);
      const std::vector <size_t> &indices = kd.find_in_range(min, max);
      for (size_t k = 0; k != indices.size(); ++k) {
        size_t j = indices[k];
        float min_distance = (radii[i] + water_radius*2) + radii[j];
        float d2 = glm::dot((pos[i] - pos[j]), (pos[i] - pos[j]));
        if (d2 * d2 <= min_distance * min_distance) {
          count++;
        }
      }

      for (size_t j = i+1; j != pos.size(); ++j) {
        float min_distance = (radii[i] + water_radius*2) + radii[j];
        float d2 = glm::dot((pos[i] - pos[j]), (pos[i] - pos[j]));
        if (d2 * d2 <= min_distance * min_distance) {
          count2++;
        }
      }
    }
    std::cout << "count = " << count << "\n";
    std::cout << "count2 = " << count2 << "\n";
    std::cout << "pos.size()*pos.size()/2 = " << pos.size()*pos.size()/2 << "\n";

    /*float spacing = 1.0f;
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
          *value++ = min_d2;
        }
      }
    }*/
  }
}
