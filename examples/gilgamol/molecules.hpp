////////////////////////////////////////////////////////////////////////////////
//
// (C) Andy Thomason 2016
//
// Example using the basic_mesh class to genererate solvent excluded
// surfaces for molecules.
//

#include <gilgamesh/mesh.hpp>
//#include <gilgamesh/distance_field.hpp>
#include <gilgamesh/decoders/pdb_decoder.hpp>
#include <gilgamesh/encoders/fbx_encoder.hpp>
#include <gilgamesh/encoders/ply_encoder.hpp>
#include <gilgamesh/shapes/sphere.hpp>
#include <gilgamesh/shapes/cylinder.hpp>
#include <gilgamesh/shapes/spline.hpp>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <ostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <thread>
#include <future>
#include <numeric>
#ifdef WIN32
  #include <io.h>
  #include <fcntl.h>
#endif

#include "utils.hpp"

class molecules {
public:
  molecules(int argc, char **argv) {
    #ifdef WIN32
      _setmode( _fileno( stdout ),  _O_BINARY );
      _setmode( _fileno( stdin ),  _O_BINARY );
    #endif
    const char *pdb_filename = nullptr;
    const char *output_path = "";
    const char *lod_text = "0";
    bool error = false;
    bool list_chains = false;
    bool use_hetatoms = false;
    const char *chains = "A-Za-z";
    const char *cmd = "";
    const char *format = "fbx";
    bool useStdio = false;

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
      } else if (!strcmp(arg, "--fbx")) {
        format = "fbx";
      } else if (!strcmp(arg, "--ply")) {
        format = "ply";
      } else if (!strcmp(arg, "--use-hetatoms")) {
        use_hetatoms = true;
      } else if (!strcmp(arg, "--list-chains")) {
        list_chains = true;
      } else if (!strcmp(arg, "-")) {
        useStdio = true;
        pdb_filename = "-";
      } else if (arg[0] == '-') {
        fprintf(stderr, "invalid argument %s\n", arg);
        error = true;
      } else if (!strcmp(arg, "se")) {
        cmd = arg;
      } else if (!strcmp(arg, "bs")) {
        cmd = arg;
      } else if (!strcmp(arg, "ca")) {
        cmd = arg;
      } else if (!strcmp(arg, "vr")) {
        cmd = arg;
      } else if (!strcmp(arg, "capsule")) {
        cmd = arg;
      } else {
        if (pdb_filename) { fprintf(stderr, "only one file will be considered\n"); error = true; }
        pdb_filename = arg;
      }
    }


    // at --lod 0, grid_spacing=1  at --lod 1, grid_spacing=0.5 etc.
    float lod = (float)atof(lod_text);
    float grid_spacing = std::pow(2.0f, -lod);

    if (pdb_filename == nullptr || error || (!cmd[0] && !list_chains)) {
      fprintf(stderr, 
        "usage:\n"
        "molecules se <options> <pdb file name> ... generate a solvent excluded mesh\n\n"
        "molecules bs <options> <pdb file name> ... generate a ball and stick mesh\n\n"
        "You can specify which chains to use with the --chains option.\n"
        "\noptions:\n"
        "--lod <n>\tLevel of detail. 0=lowest, 1=medium, 2=highest. This is logarthmic.\n"
        "--chains <n>\teg. A-E or ABDEG set of chains to use for generating FBX files. defaults to A-Z...\n"
        "--list-chains <n>\tjust list the chains in the PDB file\n"
        "--use-hetatoms\tInclude HETATM atoms\n"
        "--output-path <dir>\tdirectory to output files to\n"
        "--help <n>\tshow this text\n"
      ); return;
    }

    const char *filename = "out";

    std::vector<uint8_t> text;
    if (useStdio) {
      // 256Mb max size.
      text.resize(0x10000000);
      std::cin.read((char*)text.data(), text.size());
      size_t amount = (size_t)std::cin.gcount();
      fprintf(stderr, "amount=%d\n", (int)amount);
      text.resize(amount);
      text.shrink_to_fit();
    } else {
      std::ifstream file (pdb_filename, std::ios_base::binary);
      if (!file.eof() && !file.fail()) {
        file.seekg(0, std::ios_base::end);
        text.resize((size_t)file.tellg());

        file.seekg(0, std::ios_base::beg);
        file.read((char*)text.data(), text.size());
      }
    }
  
    gilgamesh::pdb_decoder pdb(text.data(), text.data() + text.size());

    std::string pdb_chains = pdb.chains(use_hetatoms);
  
    if (list_chains) {
      fprintf(stderr, "chains: %s\n", pdb_chains.c_str());
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

    fprintf(stderr, "using chains %s\n", expanded_chains.c_str());

    std::vector<glm::vec3> pos;
    std::vector<float> radii;
    std::vector<glm::vec4> colors;
    std::vector<std::pair<int, int> > connections;
    bool is_ca = !strcmp(cmd, "ca");
    bool is_bs = !strcmp(cmd, "bs");
    bool is_se = !strcmp(cmd, "se");
    bool is_vr = !strcmp(cmd, "vr");
    bool is_capsule = !strcmp(cmd, "capsule");

    auto atoms = pdb.atoms(expanded_chains, use_hetatoms);

    gilgamesh::color_mesh mesh;

    std::vector<glm::vec3> other_pos;
    if (is_vr) {
      // in the VR version, we colour the contact map and so need the
      // positions of the other atoms not in this chain.
      auto other_atoms = pdb.atoms(expanded_chains, true);
      for (int idx = 0; idx != other_atoms.size(); ++idx) {
        auto &p = atoms[idx];
        other_pos.push_back(glm::vec3(p.x(), p.y(), p.z()));
      }
    }

    if (is_se || is_vr) {
      for (int idx = 0; idx != atoms.size(); ++idx) {
        auto &p = atoms[idx];
        colors.push_back(p.colorByFunction());
        pos.push_back(glm::vec3(p.x(), p.y(), p.z()));
        radii.push_back(p.vanDerVaalsRadius());
      }
      generate_solvent_excluded_mesh(mesh, pos, radii, colors, other_pos, grid_spacing, is_vr);
    } else if (is_bs || is_ca) {
      for (int idx = 0; idx != atoms.size(); ++idx) {
        auto &p = atoms[idx];
        colors.push_back(p.colorByElement());
        pos.push_back(glm::vec3(p.x(), p.y(), p.z()));
        float r = p.vanDerVaalsRadius();

        if (!p.is_hetatom()) {
          // Skip atoms in alternate amino acids.
          if (p.iCode() != ' ') {
            r = 0;
          }

          // Skip 'B' alternate atoms or above
          if (p.altLoc() >= 'B') {
            r = 0;
          }

          // Ignore terminal oxygens and any Hydrogens
          if (p.atomNameIs(" OXT") || p.isHydrogen()) {
            r = 0;
          }

          if (is_ca) {
            if (p.atomNameIs(" N  ") || p.atomNameIs(" C  ") || p.atomNameIs(" O  ")) {
              r = 0;
            }
          }

          if (p.atomNameIs(" CA ") || p.atomNameIs(" N  ") || p.atomNameIs(" C  ")) {
            r *= 0.2f;
          } else {
            r *= 0.08f;
          }
        }

        radii.push_back(r);
      }

      int prevC = -1;
      char prevChainID = '?';
      for (size_t bidx = 0; bidx != atoms.size(); ) {
        // At the start of every Amino Acid, connect the atoms.
        char chainID = atoms[bidx].chainID();
        char iCode = atoms[bidx].iCode();
        size_t eidx = pdb.nextResidue(atoms, bidx);
        if (prevChainID != chainID) prevC = -1;

        // iCode is 'A' etc. for alternates.
        if (iCode == ' ') {
          prevC = pdb.addImplicitConnections(atoms, connections, bidx, eidx, prevC, is_ca);
          prevChainID = chainID;
        }
        bidx = eidx;
      }

      for (size_t i = 0; i != pos.size(); ++i) {
        bool found = false;
        for (size_t j = 0; j != connections.size() && !found; ++j) {
          found = connections[j].first == i || connections[j].second == i;
        }
        if (!found && radii[i] != 0) {
          fprintf(stderr, "atom %d %s %s not connected\n", atoms[i].serial(), atoms[i].resName().c_str(), atoms[i].atomName().c_str());
        }
      }

      generate_ball_and_stick_mesh(mesh, pos, radii, colors, connections, lod);
    } else if (is_capsule) {
      for (int idx = 0; idx != atoms.size(); ++idx) {
        auto &p = atoms[idx];
        pos.push_back(glm::vec3(p.x(), p.y(), p.z()));
      }
      generate_capsule_mesh(mesh, pos);
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

    const char *out_filename = fmt("%s_%s_%s_%s.%s", stem.c_str(), expanded_chains.c_str(), cmd, lod_text, format);
    if (useStdio) {
      out_filename = "-";
    }
    fprintf(stderr, "writing %s (%d vertices)\n", out_filename, int(mesh.vertices().size()));

    if (format[0] == 'p') {
      gilgamesh::ply_encoder encoder;
      encoder.saveMesh(mesh, out_filename, false, "pnuc");
    } else {
      gilgamesh::fbx_encoder encoder;
      encoder.saveMesh(mesh, out_filename);
    }
  }

private:
  void generate_ball_and_stick_mesh(gilgamesh::color_mesh &mesh, std::vector<glm::vec3> &pos, std::vector<float> &radii, std::vector<glm::vec4> &colors, std::vector<std::pair<int, int> > &connections,float lod) {
    glm::mat4 mat;
    for (size_t i = 0; i != pos.size(); ++i) {
      if (radii[i] != 0) {
        mat[3].x = pos[i].x; mat[3].y = pos[i].y; mat[3].z = pos[i].z;
        gilgamesh::sphere s(radii[i]);
        int segments = int(5 + lod * 4);
        s.build(mesh, mat, colors[i], segments);
      }
    }

    for (auto &c : connections) {
      //fprintf(stderr, "%d %d / %d\n", c.first, c.second, (int)pos.size());
      glm::vec3 pos0 = pos[c.first];
      glm::vec3 pos1 = pos[c.second];
      glm::vec4 c0 = colors[c.first];
      glm::vec4 c1 = colors[c.second];
      glm::vec3 up = glm::vec3(0, 1, 0);
      glm::vec3 y = glm::normalize(pos1 - pos0);
      glm::vec3 x = glm::normalize(glm::cross(y, up));
      glm::vec3 z = glm::cross(x, y);
      float len = glm::length(pos1 - pos0);
      if (len < 5) {
        mat[0] = glm::vec4(x, 0);
        mat[1] = glm::vec4(y, 0);
        mat[2] = glm::vec4(z, 0);
        /*if (len >= 3.0 || len <= 0.001f) {
          fprintf(stderr, "%d %d %f %f %f  %f %f %f %f\n", c.first, c.second, pos0.x, pos0.y, pos0.z, mat[3].x, mat[3].y, mat[3].z, len);
        }*/
        float r = std::min(radii[c.first], radii[c.second]) * 0.5f;
        int segments = 5;
        if (c0 == c1) {
          gilgamesh::cylinder cyl(r, len);
          mat[3] = glm::vec4((pos0 + pos1) * 0.5f, 1);
          cyl.build(mesh, mat, colors[c.first], 1, segments, gilgamesh::cylinder::body);
        } else {
          gilgamesh::cylinder cyl(r, len * 0.5f);
          mat[3] = glm::vec4(pos0 * 0.75f + pos1 * 0.25f, 1);
          cyl.build(mesh, mat, colors[c.first], 1, segments, gilgamesh::cylinder::body);
          mat[3] = glm::vec4(pos0 * 0.25f + pos1 * 0.75f, 1);
          cyl.build(mesh, mat, colors[c.second], 1, segments, gilgamesh::cylinder::body);
        }
      }
    }
  }

  void generate_solvent_excluded_mesh(gilgamesh::color_mesh &mesh, std::vector<glm::vec3> &pos, std::vector<float> &radii, std::vector<glm::vec4> &colors, std::vector<glm::vec3> &other_pos, float grid_spacing, bool is_vr) {
    struct colored_atom {
      glm::vec4 color;
      glm::vec3 pos;
      float radius;
    };
    std::vector<colored_atom> colored_atoms;

    glm::vec4 white(1, 1, 1, 1);
    glm::vec4 red(1, 0, 0, 1);
    glm::vec4 blue(0, 0, 1, 1);
    for (size_t i = 0; i != colors.size(); ++i) {
      if (colors[i] == red || colors[i] == blue) {
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
    fprintf(stderr, "%d x %d x %d\n", xdim, ydim, zdim);

    auto idx = [xdim, ydim](int x, int y, int z) {
      return ((z * (ydim+1)) + y) * (xdim+1) + x;
    };

    fprintf(stderr, "building solvent acessible mesh by inflating the atoms\n");

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

    fprintf(stderr, "building solvent excluded mesh by deflating the acessible mesh\n");
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

      bool qq = false;
      float tot = 1;
      for (size_t i = 0; i != colored_atoms.size(); ++i) {
        glm::vec3 &pos = colored_atoms[i].pos;
        float r = colored_atoms[i].radius;
        float r2 = r * r;
        float d2 = glm::dot(xyz - pos, xyz - pos);
        float e0 = r2 * 16;
        if (d2 < e0) {
          qq = true;
          float e1 = r2 * 4;
          float v = std::min(std::max((d2 - e0)/(e1 - e0), 0.0f), 1.0f);
          float weight = v * v * (3 - 2 * v);
          color += colored_atoms[i].color * weight;
          tot += weight;
          //fprintf(stderr, "%f %f %f %f\n", colored_atoms[i].color.x, colored_atoms[i].color.y, colored_atoms[i].color.z, weight);
          //fprintf(stderr, "%f %f %f\n", color.x, color.y, color.z);
        }
      }

      color /= tot;

      color.x = std::min(std::max(color.x, 0.0f), 1.0f);
      color.y = std::min(std::max(color.y, 0.0f), 1.0f);
      color.z = std::min(std::max(color.z, 0.0f), 1.0f);

      return gilgamesh::color_mesh::vertex_t(xyz, normal, uv, color);
    };

    auto vr_egen = [&excluded, &other_pos, grid_spacing, min, idx](float x, float y, float z) {
      glm::vec3 xyz(x * grid_spacing + min.x, y * grid_spacing + min.y, z * grid_spacing + min.z);
      glm::vec3 normal(1, 0, 0);
      glm::vec2 uv(0, 0);
      glm::vec4 color = glm::vec4(1, 1, 1, 1);

      auto low = std::lower_bound(other_pos.begin(), other_pos.end(), x - 4, [](const glm::vec3 &a, float value) { return a.x < value; });
      auto high = std::upper_bound(other_pos.begin(), other_pos.end(), x + 4, [](float value, const glm::vec3 &a) { return value < a.x; });
      fprintf(stderr, "%f %f %f %d %d\n", x, y, z, (int)(low - other_pos.begin()), (int)(high - other_pos.begin()));

      glm::vec3 pos(x, y, z);
      for (auto p = low; p != high; ++p) {
        glm::vec3 d = *p - pos;
        float d2 = glm::dot(d, d);
        fprintf(stderr, "  %f %f %f %f\n", p->x, p->y, p->z, d2);
        if (d2 < 16) {
          color = glm::vec4(0.5f, 0.5f, 0.5f, 1);
          fprintf(stderr, "hooray!\n");
        }
      }

      return gilgamesh::color_mesh::vertex_t(xyz, normal, uv, color);
    };

    // construct the excluded mesh using marching cubes
    if (is_vr) {
      std::sort(other_pos.begin(), other_pos.end(), [](const glm::vec3 &a, const glm::vec3 &b) { return a.x < b.x; });
      mesh = gilgamesh::color_mesh(xdim, ydim, zdim, efn, vr_egen);
    } else {
      mesh = gilgamesh::color_mesh(xdim, ydim, zdim, efn, egen);
    }
  }

  void generate_capsule_mesh(gilgamesh::mesh &mesh, const std::vector<glm::vec3> &pos) {
    /*gilgamesh::Spline spline(pos, gilgamesh::SplineType::CatmullRom);

    gilgamesh::simple_mesh splinePoints;
    spline.build(splinePoints);

    gligamesh::Circle circle(1.0f);
    gilgamesh::simple_mesh circlePoints;
    circle.build(circlePoints);

    gilgamesh::Loft loft(splinePoints.pos(), splinePoints.uvs(), circlePoints.pos(), circlePoints.uvs());
    circle.build(mesh);*/
  }
};

