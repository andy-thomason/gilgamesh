////////////////////////////////////////////////////////////////////////////////
//
// (C) Andy Thomason 2016
//
// Example using the basic_mesh class to genererate solvent excluded
// surfaces for molecules.
//

#include <gilgamesh/mesh.hpp>
#include <gilgamesh/decoders/pdb_decoder.hpp>
#include <gilgamesh/encoders/fbx_encoder.hpp>
#include <gilgamesh/shapes/sphere.hpp>
#include <gilgamesh/shapes/cylinder.hpp>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <ostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <thread>
#include <future>
#include <numeric>

#include "utils.hpp"

class molecules {
public:
  molecules(int argc, char **argv) {
    const char *pdb_filename = nullptr;
    const char *output_path = "";
    const char *lod_text = "0";
    bool error = false;
    bool list_chains = false;
    const char *chains = "A-Z";
    const char *cmd = "";

    for (int i = 1; i != argc; ++i) {
      const char *arg = argv[i];
      if (!strcmp(arg, "--lod") && i < argc-1) {
        lod_text = argv[++i];
      } else if (!strcmp(arg, "--chains") && i < argc-1) {
        chains = argv[++i];
      } else if (!strcmp(arg, "-o") && i < argc-1) {
        output_path = argv[++i];
      } else if (!strcmp(arg, "--help")) {
        error = true;
      } else if (!strcmp(arg, "--list-chains")) {
        list_chains = true;
      } else if (arg[0] == '-') {
        printf("invalid argument %s\n", arg);
        error = true;
      } else if (!strcmp(arg, "se")) {
        cmd = arg;
      } else if (!strcmp(arg, "bs")) {
        cmd = arg;
      } else if (!strcmp(arg, "ca")) {
        cmd = arg;
      } else {
        if (pdb_filename) { printf("only one file will be considered\n"); error = true; }
        pdb_filename = arg;
      }
    }

    // at --lod 0, grid_spacing=1  at --lod 1, grid_spacing=0.5 etc.
    float grid_spacing = std::pow(2.0f, -float(atof(lod_text)));

    if (pdb_filename == nullptr || error || !cmd[0]) {
      printf(
        "usage:\n"
        "molecules se <options> <pdb file name> ... generate a solvent excluded mesh\n\n"
        "molecules bs <options> <pdb file name> ... generate a ball and stick mesh\n\n"
        "You can specify which chains to use with the --chains option.\n"
        "\noptions:\n"
        "--lod <n>\tLevel of detail. 0=lowest, 1=medium, 2=highest. This is logarthmic.\n"
        "--chains <n>\teg. A-E or ABDEG set of chains to use for generating FBX files. defaults to A-Z...\n"
        "--list-chains <n>\tjust list the chains in the PDB file\n"
        "--output-path <dir>\tdirectory to output files to\n"
        "--help <n>\tshow this text\n"
      ); return;
    }

    const char *filename = "out";

    //std::ifstream file(CMAKE_SOURCE "/examples/data/2PTC.pdb", std::ios_base::binary);
    std::ifstream file(pdb_filename, std::ios_base::binary);
    std::vector<uint8_t> text;
    if (!file.eof() && !file.fail()) {
      file.seekg(0, std::ios_base::end);
      text.resize((size_t)file.tellg());

      file.seekg(0, std::ios_base::beg);
      file.read((char*)text.data(), text.size());
    }
  
    gilgamesh::pdb_decoder pdb(text.data(), text.data() + text.size());

    gilgamesh::fbx_encoder encoder;
    std::string pdb_chains = pdb.chains();
  
    if (list_chains) {
      printf("chains: %s\n", pdb_chains.c_str());
      return;
    }

    std::string expanded_chains;
    for (const char *p = chains; *p; ++p) {
      if (p[1] == '-' && p[2] && p[2] >= p[1]) {
        for (char i = p[0]; i <= p[2]; ++i) {
          if (pdb_chains.find(i) != std::string::npos) {
            expanded_chains.push_back(i);
          }
        }
        p += 2;
      } else {
        expanded_chains.push_back(*p);
      }
    }

    printf("chains %s\n", expanded_chains.c_str());

    std::vector<glm::vec3> pos;
    std::vector<float> radii;
    std::vector<glm::vec4> colors;
    std::vector<std::pair<int, int> > connections;
    bool is_ca = !strcmp(cmd, "ca");
    bool is_bs = !strcmp(cmd, "bs");
    bool is_se = !strcmp(cmd, "se");

    auto atoms = pdb.atoms();
    int prevC = -1;
    char prevChainID = '?';
    int prev_resSeq = -1;
    for (int idx = 0; idx != atoms.size(); ++idx) {
      auto &p = atoms[idx];
      char chainID = p.chainID();
      // if the chain is in the set specified on the command line (eg. ACBD)
      if (expanded_chains.find(chainID) != std::string::npos) {
        // if we are not in CA only mode or we are in the set of CA atoms (C, N, CA)
        if (!is_ca || p.atomNameIs(" C  ") || p.atomNameIs(" N  ") || p.atomNameIs(" CA ")) {
          if (is_bs || is_ca) {
            colors.push_back(p.colorByElement());
            if (p.resSeq() != prev_resSeq) {
              int resSeq = p.resSeq(), j = 0;
              prev_resSeq = resSeq;
              for (j = idx; j != atoms.size(); ++j) {
                if (atoms[j].resSeq() != resSeq) {
                  break;
                }
              }
              const auto *b = atoms.data() + idx;
              const auto *e = atoms.data() + j;

              // At the start of every Amino Acid, connect the atoms.
              int N_idx = int(pos.size());
              int C_idx = gilgamesh::pdb_decoder::addImplicitConnections(connections, b, e, N_idx, is_ca);
              if (prevC != -1 && prevChainID == chainID) {
                connections.emplace_back(prevC, N_idx);
              }
              prevC = C_idx;
              prevChainID = chainID;
            }
          } else {
            colors.push_back(p.colorByFunction());
          }
          pos.push_back(glm::vec3(p.x(), p.y(), p.z()));
          radii.push_back(p.vanDerVaalsRadius());
        }
      }
    }

    gilgamesh::color_mesh mesh;

    if (is_se) {
      generate_solvent_excluded_mesh(mesh, pos, radii, colors, grid_spacing);
    }

    if (is_bs || is_ca) {
      generate_ball_and_stick_mesh(mesh, pos, radii, colors, connections);
    }

    const char *last_slash = pdb_filename;
    const char *last_dot = pdb_filename + strlen(pdb_filename);
    for (const char *p = pdb_filename; *p; ++p) {
      if (*p == '/' || *p == '\\') last_slash = p + 1;
    }
    for (const char *p = last_slash; *p; ++p) {
      if (*p == '.') last_dot = p;
    }

    // build normals
    mesh.reindex(true);

    std::string stem;
    stem.assign(last_slash, last_dot);

    const char *out_filename = fmt("%s_%s_%s_%s.fbx", stem.c_str(), expanded_chains.c_str(), cmd, lod_text);
    printf("writing %s (%d vertices)\n", out_filename, int(mesh.vertices().size()));
    std::ofstream eof(out_filename, std::ios::binary);
    std::vector<uint8_t> bytes = encoder.saveMesh(mesh);
    eof.write((char*)bytes.data(), bytes.size());
  }

private:
  void generate_ball_and_stick_mesh(gilgamesh::color_mesh &mesh, std::vector<glm::vec3> &pos, std::vector<float> &radii, std::vector<glm::vec4> &colors, std::vector<std::pair<int, int> > &connections) {
    glm::mat4 mat;
    for (size_t i = 0; i != pos.size(); ++i) {
      mat[3].x = pos[i].x; mat[3].y = pos[i].y; mat[3].z = pos[i].z;
      gilgamesh::sphere s(0.40f);
      s.build(mesh, mat, colors[i], 5);
    }

    for (auto &c : connections) {
      glm::vec3 pos0 = pos[c.first];
      glm::vec3 pos1 = pos[c.second];
      glm::vec4 c0 = colors[c.first];
      glm::vec4 c1 = colors[c.second];
      glm::vec3 up = glm::vec3(0, 1, 0);
      glm::vec3 y = glm::normalize(pos1 - pos0);
      glm::vec3 x = glm::normalize(glm::cross(y, up));
      glm::vec3 z = glm::cross(x, y);
      float len = glm::length(pos1 - pos0);
      if (len < 3.0f) {
        mat[0] = glm::vec4(x, 0);
        mat[1] = glm::vec4(y, 0);
        mat[2] = glm::vec4(z, 0);
        /*if (len >= 3.0 || len <= 0.001f) {
          printf("%d %d %f %f %f  %f %f %f %f\n", c.first, c.second, pos0.x, pos0.y, pos0.z, mat[3].x, mat[3].y, mat[3].z, len);
        }*/
        if (c0 == c1) {
          gilgamesh::cylinder cyl(0.25f, len);
          mat[3] = glm::vec4((pos0 + pos1) * 0.5f, 1);
          cyl.build(mesh, mat, colors[c.first], 1, 8, gilgamesh::cylinder::body);
        } else {
          gilgamesh::cylinder cyl(0.25f, len * 0.5f);
          mat[3] = glm::vec4(pos0 * 0.75f + pos1 * 0.25f, 1);
          cyl.build(mesh, mat, colors[c.first], 1, 8, gilgamesh::cylinder::body);
          mat[3] = glm::vec4(pos0 * 0.25f + pos1 * 0.75f, 1);
          cyl.build(mesh, mat, colors[c.second], 1, 8, gilgamesh::cylinder::body);
        }
      }
    }
  }

  void generate_solvent_excluded_mesh(gilgamesh::color_mesh &mesh, std::vector<glm::vec3> &pos, std::vector<float> &radii, std::vector<glm::vec4> &colors, float grid_spacing) {
    struct colored_atom {
      glm::vec4 color;
      glm::vec3 pos;
      float radius;
    };
    std::vector<colored_atom> colored_atoms;

    glm::vec4 white(1, 1, 1, 1);
    for (size_t i = 0; i != colors.size(); ++i) {
      if (colors[i] != white) {
        colored_atom a = { colors[i], pos[i], radii[i] };
        colored_atoms.push_back(a);
      }
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

    printf("building solvent acessible mesh by inflating the atoms\n");

    // build a distance field for all the atoms.
    std::vector<float> accessible((xdim+1)*(ydim+1)*(zdim+1));
    par_for(0, zdim+1, [&](int z) {
      float zpos = z * grid_spacing + min.z;
      for (int y = 0; y != ydim+1; ++y) {
        float ypos = y * grid_spacing + min.y;
        for (int x = 0; x != xdim+1; ++x) {
          glm::vec3 xyz(x * grid_spacing + min.x, ypos, zpos);
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
      return gilgamesh::pos_mesh::vertex_t(xyz);
    };

    // construct the accessible mesh using marching cubes
    gilgamesh::pos_mesh amesh(xdim, ydim, zdim, fn, gen);

    std::vector<glm::vec3> zsorter;

    auto &avertices = amesh.vertices();
    for (auto &v : avertices) {
      zsorter.push_back(v.pos());
    }

    auto cmpz = [](const glm::vec3 &a, const glm::vec3 &b) { return a.z < b.z; };
    auto cmpy = [](const glm::vec3 &a, const glm::vec3 &b) { return a.y < b.y; };
    std::sort(zsorter.begin(), zsorter.end(), cmpz);

    printf("building solvent excluded mesh by deflating the acessible mesh\n");
    std::vector<float> excluded((xdim+1)*(ydim+1)*(zdim+1));
    float outside_value = -(water_radius * water_radius);
    par_for(0, zdim+1, [&](int z) {
      std::vector<glm::vec3> ysorter;

      // search only a band of z values in zpos +/- water_radius
      float zpos = z * grid_spacing + min.z;
      auto p = std::lower_bound( zsorter.begin(), zsorter.end(), glm::vec3(0, 0, zpos - (water_radius + grid_spacing)), cmpz);
      auto q = std::upper_bound( zsorter.begin(), zsorter.end(), glm::vec3(0, 0, zpos + (water_radius + grid_spacing)), cmpz);

      for (int y = 0; y != ydim+1; ++y) {
        float ypos = y * grid_spacing + min.y;

        // filter by y position
        ysorter.clear();
        for (auto r = p; r != q; ++r) {
          if (std::abs(r->y - ypos) <= water_radius + grid_spacing) {
            ysorter.push_back(*r);
          }
        }

        for (int x = 0; x != xdim+1; ++x) {
          float value = 1e37f;
          // only if we are inside the acessible mesh...
          if (accessible[idx(x, y, z)] < 0) {
            glm::vec3 xyz(x * grid_spacing + min.x, ypos, zpos);
            // find the closest point on the accessible mesh to xyz.
            for (auto &r : ysorter) {
              float d2 = glm::dot(xyz - r, xyz - r);
              value = std::min(value, d2);
            }
            if (value == 1e37f) {
              value = outside_value;
            } else {
              value -= (water_radius) * (water_radius);
            }
          } else {
            value = outside_value;
          }
          excluded[idx(x, y, z)] = value;
        }
      }
    });

    auto efn = [&excluded, idx](int x, int y, int z) {
      return excluded[idx(x, y, z)];
    };

    auto egen = [&excluded, &colored_atoms, grid_spacing, min, idx](float x, float y, float z) {
      glm::vec3 xyz(x * grid_spacing + min.x, y * grid_spacing + min.y, z * grid_spacing + min.z);
      glm::vec3 normal(1, 0, 0);
      glm::vec2 uv(0, 0);
      glm::vec4 color = glm::vec4(1, 1, 1, 1);
      int num_influences = 1;
      for (size_t i = 0; i != colored_atoms.size(); ++i) {
        glm::vec3 &pos = colored_atoms[i].pos;
        float r = colored_atoms[i].radius;
        float d2 = glm::dot(xyz - pos, xyz - pos);
        float weight = (r*r*4 - d2) * 1.0f;
        if (weight > 0) {
          weight = std::max(0.0f, std::min(weight, 1.0f));
          color += colored_atoms[i].color * weight;
          num_influences++;
        }
      }
      color.w = 1;
      color.x *= (1.0f/num_influences);
      color.y *= (1.0f/num_influences);
      color.z *= (1.0f/num_influences);
      return gilgamesh::color_mesh::vertex_t(xyz, normal, uv, color);
    };

    // construct the excluded mesh using marching cubes
    gilgamesh::color_mesh emesh(xdim, ydim, zdim, efn, egen);

    // use move operator to shallow copy the mesh.
    mesh = std::move(emesh);
  }
};

