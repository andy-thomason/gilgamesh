////////////////////////////////////////////////////////////////////////////////
//
// (C) Andy Thomason 2016
//
// Fbx file decoder
// 

#ifndef VKU_fbx_decoder_INCLUDED
#define VKU_fbx_decoder_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
//#include <filesystem>
#include <vector>

#include <meshutils/scene.hpp>
#include <minizip/deflate_decoder.hpp>
#include <glm/glm.hpp>

// see https://code.blender.org/2013/08/fbx-binary-file-format-specification/
// and https://banexdevblog.wordpress.com/2014/06/23/a-quick-tutorial-about-the-fbx-ascii-format/

namespace meshutils {

  class fbx_decoder {
    enum { debug = 1 };

    static inline std::uint8_t u1(const char *p) {
      const unsigned char *q = (const unsigned char *)p;
      std::uint32_t res = q[0];
      return res;
    }

    static inline std::uint16_t u2(const char *p) {
      const unsigned char *q = (const unsigned char *)p;
      std::uint32_t res = q[0] + q[1] * 0x100;
      return res;
    }

    static inline std::uint32_t u4(const char *p) {
      const unsigned char *q = (const unsigned char *)p;
      std::uint32_t res = q[0] + q[1] * 0x100 + q[2] * 0x10000 + q[3] * 0x1000000;
      return res;
    }

    static inline std::uint64_t u8(const char *p) {
      return ((std::uint64_t)u4(p+4) << 32) | u4(p);
    }

    class prop;

    class props {
    public:
      props(const char *begin=nullptr, size_t offset=0) : begin_(begin), offset(offset) {}
      prop begin() const { return prop(begin_, offset + 13 + len()); }
      prop end() const { return prop(begin_, offset + 13 + property_list_len() + len()); }

    private:
      size_t end_offset() const { return u4(begin_ + offset); }
      size_t num_properties() const { return u4(begin_ + offset + 4); }
      size_t property_list_len() const { return u4(begin_ + offset + 8); }
      std::uint8_t len() const { return u1(begin_ + offset + 12); }

      size_t offset;
      const char *begin_;
    };

    // http://code.blender.org/2013/08/fbx-binary-file-format-specification/
    class node {
    public:
      node(const char *begin=nullptr, size_t offset=0) : begin_(begin), offset_(offset) {}
      bool operator !=(node &rhs) { return offset_ != rhs.offset_; }
      node &operator++() { offset_ = end_offset(); return *this; }
      std::string name() const { return std::string(begin_ + offset_ + 13, begin_ + offset_ + 13 + len()); }
      node begin() const { size_t new_offset = offset_ + 13 + property_list_len() + len(), end = end_offset(); return node(begin_, new_offset == end ? end-13 : new_offset); }
      node end() const { return node(begin_, end_offset() - 13); }
      node &operator*() { return *this; }
      fbx_decoder::props get_props() { return fbx_decoder::props(begin_, offset_); }

      size_t offset() const { return offset_; }
      size_t end_offset() const { return u4(begin_ + offset_); }
      size_t num_properties() const { return u4(begin_ + offset_ + 4); }
      size_t property_list_len() const { return u4(begin_ + offset_ + 8); }
      std::uint8_t len() const { return u1(begin_ + offset_ + 12); }

    //private:
      size_t offset_;
      const char *begin_;
    };

    class prop {
    public:
      prop(const char *begin=nullptr, size_t offset=0) : begin_(begin), offset(offset) {}
      bool operator !=(prop &rhs) { return offset != rhs.offset; }
      prop &operator*() { return *this; }
      char kind() const { return begin_[offset]; }

      prop &operator++() {
        const char *p = begin_ + offset;
        size_t al, enc, cl;
        switch(*p++) {
            case 'Y': p += 2; break;
            case 'C': p += 1; break;
            case 'I': p += 4; break;
            case 'F': p += 4; break;
            case 'D': p += 8; break;
            case 'L': p += 8; break;
            case 'f': al = u4(p); enc = u4(p+4); cl = u4(p+8); p += 12; p += !enc ? al * 4 : cl; break;
            case 'd': al = u4(p); enc = u4(p+4); cl = u4(p+8); p += 12; p += !enc ? al * 8 : cl; break;
            case 'l': al = u4(p); enc = u4(p+4); cl = u4(p+8); p += 12; p += !enc ? al * 8 : cl; break;
            case 'i': al = u4(p); enc = u4(p+4); cl = u4(p+8); p += 12; p += !enc ? al * 4 : cl; break;
            case 'b': al = u4(p); enc = u4(p+4); cl = u4(p+8); p += 12; p += !enc ? al * 1 : cl; break;
            case 'S': al = u4(p); p += 4 + al; break;
            case 'R': al = u4(p); p += 4 + al; break;
            default: throw std::runtime_error("bad fbx property"); break;
        }
        offset = p - begin_;
        return *this;
      }

      char *ilist(char *d, char *e, const char *p, size_t al, size_t enc, size_t cl, size_t elem_size) {
        d += snprintf(d, e-d, "%d, %d, %d, {", (int)al, (int)enc, (int)cl);
        size_t size = enc ? cl : al * elem_size;
        while (d < e && size--) {
          int c = *p++ & 0xff;
          d += snprintf(d, e-d, "0x%02x,", c);
        }
        d += snprintf(d, e-d, "}");
        return d;
      }

      operator std::string() {
        const char *p = begin_ + offset;
        size_t al, enc, cl;
        static char tmp[65536];
        char *d = tmp, *e = tmp + sizeof(tmp) - 10;
        int fv; std::uint64_t dv;
        switch(*p++) {
          case 'Y': snprintf(tmp, sizeof(tmp), "%d", (short)u2(p)); break;
          case 'C': snprintf(tmp, sizeof(tmp), *p ? "true" : "false"); break;
          case 'I': snprintf(tmp, sizeof(tmp), "%d", (std::int32_t)u4(p)); break;
          case 'F': fv = u4(p); snprintf(tmp, sizeof(tmp), "%8f", (float&)(fv)); break;
          case 'D': dv = u8(p); snprintf(tmp, sizeof(tmp), "%10f", (double&)(dv)); break;
          case 'L': snprintf(tmp, sizeof(tmp), "%lld", (long long)u8(p)); break;
          case 'f': al = u4(p); enc = u4(p+4); cl = u4(p+8); d = ilist(d, e, p + 12, al, enc, cl, 4); break;
          case 'd': al = u4(p); enc = u4(p+4); cl = u4(p+8); d = ilist(d, e, p + 12, al, enc, cl, 8); break;
          case 'l': al = u4(p); enc = u4(p+4); cl = u4(p+8); d = ilist(d, e, p + 12, al, enc, cl, 8); break;
          case 'i': al = u4(p); enc = u4(p+4); cl = u4(p+8); d = ilist(d, e, p + 12, al, enc, cl, 4); break;
          case 'b': al = u4(p); enc = u4(p+4); cl = u4(p+8); d = ilist(d, e, p + 12, al, enc, cl, 4); break;
          case 'S': {
            size_t size = al = u4(p);
            bool has_null = false;
            *d++ = '"';
            p += 4;
            while (d < e && al--) {
              int c = *p++ & 0xff;
              if (c < ' ' || c == '\\' || c == '"' || c >= 0x7f) {
                if (c == 0) has_null = true;
                d += snprintf(d, e-d, "\\x%02x", c);
              } else {
                *d++ = c;
              }
            }
            *d++ = '"';
            if (has_null) {
              d += snprintf(d, e-d, ", %d", (int)size);
            }
            *d = 0;
          } break;
          case 'R': al = u4(p); d = ilist(d, e, p + 4, al, 0, al, 1); break;
          default: throw std::runtime_error("bad fbx property"); break;
        }
        return tmp;
      }

      bool getString(std::string &result) {
        if (kind() == 'S') {
          const char *p = begin_ + offset + 1;
          size_t al = u4(p);
          result.assign(p+4, p+4+al);
          return true;
        }
        return false;
      }

      template <class Type, char Kind>
      bool getArray(std::vector<Type> &result, const minizip::deflate_decoder &decoder) const {
        Type *begin = nullptr;
        Type *end = nullptr;
        if (kind() == Kind) {
          const char *p = begin_ + offset + 1;
          size_t al = u4(p);
          size_t enc = u4(p+4);
          size_t cl = u4(p+8);
          p += 12;
          result.resize(al);
          uint8_t *dest = (uint8_t *)result.data();
          uint8_t *dest_max = (uint8_t *)(result.data() + result.size());
          if (enc) {
            const uint8_t *src = (const uint8_t *)p;
            const uint8_t *src_max = (const uint8_t *)p + cl;
            // bytes 0 and 1 are the ZLIB code.
            // see http://stackoverflow.com/questions/9050260/what-does-a-zlib-header-look-like
            if ((src[0] & 0x0f) == 0x08) {
              return decoder.decode(dest, dest_max, src+2, src_max);
            }
          } else {
            memcpy(dest, p, al*sizeof(Type));
            return true;
          }
        }
        return false;
      }
    private:
      size_t offset;
      const char *begin_;
    };
  public:
    /*fbx_decoder(const std::string &filename) {
      std::ifstream f(filename, std::ios_base::binary);
      if (f.good()) {
        f.seekg(0, std::ios_base::end);
        bytes.resize((size_t)f.tellg());
        f.seekg(0, std::ios_base::beg);
        f.read((char*)bytes.data(), bytes.size());
        init((char*)bytes.data(), (char*)bytes.data() + f.gcount());
      }
    }*/

    fbx_decoder(const char *begin, const char *end) { init(begin, end); }

    node begin() const { return node(begin_, 27); }
    node end() const { return node(begin_, end_offset); }

    enum class Mapping {
      Invalid,
      ByPolygon,
      ByPolygonVertex,
      ByVertex,
      ByEdge,
      AllSame,
    };

    enum class Ref {
      Invalid,
      Direct,
      IndexToDirect,
    };

    static Mapping decodeMapping(const std::string &name) {
      Mapping result = Mapping::Invalid;
      if (name == "ByPolygon") result = Mapping::ByPolygon;
      else if (name == "ByPolygon") result = Mapping::ByPolygon;
      else if (name == "ByPolygonVertex") result = Mapping::ByPolygonVertex;
      else if (name == "ByVertex") result = Mapping::ByVertex;
      else if (name == "ByVertice") result = Mapping::ByVertex;
      else if (name == "ByEdge") result = Mapping::ByEdge;
      else if (name == "AllSame") result = Mapping::AllSame;
      return result;
    };

    static Ref decodeRef(const std::string &name) {
      Ref result = Ref::Invalid;
      if (name == "Direct") result = Ref::Direct;
      else if (name == "IndexToDirect") result = Ref::IndexToDirect;
      else if (name == "Index") result = Ref::IndexToDirect;
      return result;
    };

    template<class MeshType>
    bool loadScene(meshutils::scene &scene) {
      minizip::deflate_decoder decoder;

      std::vector<double> fbxVertices;
      std::vector<double> fbxNormals;
      std::vector<double> fbxUVs;
      std::vector<int32_t> fbxUVIndices;
      std::vector<int32_t> fbxNormalIndices;
      std::vector<int32_t> fbxIndices;
      std::string fbxNormalMapping;
      std::string fbxUVMapping;
      std::string fbxNormalRef;
      std::string fbxUVRef;

      for (auto section : *this) {
        if (section.name() == "Objects") {
          for (auto obj : section) {
            if (obj.name() == "Geometry") {
              for (auto comp : obj) {
                auto vp = comp.get_props().begin();
                if (debug) printf("%s %c\n", comp.name().c_str(), vp.kind());
                if (comp.name() == "Vertices") {
                  vp.getArray<double, 'd'>(fbxVertices, decoder);
                } else if (comp.name() == "LayerElementNormal") {
                  for (auto sub : comp) {
                    auto vp = sub.get_props().begin();
                    if (debug) printf("  %s %c\n", sub.name().c_str(), vp.kind());
                    if (sub.name() == "MappingInformationType") {
                      vp.getString(fbxNormalMapping);
                    } else if (sub.name() == "ReferenceInformationType") {
                      vp.getString(fbxNormalRef);
                    } else if (sub.name() == "NormalIndex") {
                      vp.getArray<int32_t, 'i'>(fbxNormalIndices, decoder);
                    } else if (sub.name() == "Normals") {
                      vp.getArray<double, 'd'>(fbxNormals, decoder);
                    }
                  }
                } else if (comp.name() == "LayerElementUV") {
                  for (auto sub : comp) {
                    auto vp = sub.get_props().begin();
                    if (debug) printf("  %s %c\n", sub.name().c_str(), vp.kind());
                    if (sub.name() == "MappingInformationType") {
                      vp.getString(fbxUVMapping);
                    } else if (sub.name() == "ReferenceInformationType") {
                      vp.getString(fbxUVRef);
                    } else if (sub.name() == "UVIndex") {
                      vp.getArray<int32_t, 'i'>(fbxUVIndices, decoder);
                    } else if (sub.name() == "UV") {
                      vp.getArray<double, 'd'>(fbxUVs, decoder);
                    }
                  }
                } else if (comp.name() == "PolygonVertexIndex") {
                  vp.getArray<int32_t, 'i'>(fbxIndices, decoder);
                }
              }

              auto normalMapping = fbx_decoder::decodeMapping(fbxNormalMapping);
              auto uvMapping = fbx_decoder::decodeMapping(fbxUVMapping);
              auto normalRef = fbx_decoder::decodeRef(fbxNormalRef);
              auto uvRef = fbx_decoder::decodeRef(fbxUVRef);

              if (fbxNormals.empty()) {
                fbxNormals.resize(3);
              }
              if (fbxUVs.empty()) {
                fbxUVs.resize(2);
              }

              // https://banexdevblog.wordpress.com/2014/06/23/a-quick-tutorial-about-the-fbx-ascii-format/
              if (debug) printf("%s %s\n", fbxNormalMapping.c_str(), fbxUVMapping.c_str());
              if (debug) printf("%s %s\n", fbxNormalRef.c_str(), fbxUVRef.c_str());
              if (debug) printf("%d vertices %d indices %d normals %d uvs %d uvindices\n", (int)fbxVertices.size(), (int)fbxIndices.size(), (int)fbxNormals.size(), (int)fbxUVs.size(), (int)fbxUVIndices.size());

              std::vector<glm::vec3> pos;
              std::vector<glm::vec3> normal;
              std::vector<glm::vec2> uv;
              std::vector<glm::vec4> color;
              std::vector<int> material;

              // map the fbx data to real vertices
              size_t pi = 0;
              for (size_t i = 0; i != fbxIndices.size(); ++i) {
                size_t ni = normalRef == fbx_decoder::Ref::IndexToDirect ? fbxNormalIndices[i] : i;
                size_t uvi = uvRef == fbx_decoder::Ref::IndexToDirect ? fbxUVIndices[i] : i;
                int32_t vi = fbxIndices[i];
                if (vi < 0) vi = -1 - vi;

                glm::vec3 vpos(fbxVertices[vi*3+0], fbxVertices[vi*3+1], fbxVertices[vi*3+2]);
                glm::vec3 vnormal(1, 0, 0);
                glm::vec2 vuv(0, 0);
                glm::vec4 vcolor(1, 1, 1, 1);

                size_t nj = map(normalMapping, pi, ni, vi);
                size_t uvj = map(uvMapping, pi, uvi, vi);

                vnormal = glm::vec3(fbxNormals[nj*3+0], fbxNormals[nj*3+1], fbxNormals[nj*3+2]);
                vuv = glm::vec2(fbxUVs[uvj*2+0], fbxUVs[uvj*2+1]);

                pos.push_back(vpos);
                normal.push_back(vnormal);
                uv.push_back(vuv);
                color.push_back(vcolor);

                pi += fbxIndices[i] < 0;
              }


              // map the fbx data to real indices
              // todo: add a function to re-index
              std::vector<uint32_t> indices;
              for (size_t i = 0, j = 0; i != fbxIndices.size(); ++i) {
                if (fbxIndices[i] < 0) {
                  for (size_t k = j+2; k <= i; ++k) {
                    indices.push_back((uint32_t)j);
                    indices.push_back((uint32_t)k-1);
                    indices.push_back((uint32_t)k);
                  }
                  j = i + 1;
                }
              }

              MeshType *mesh = new MeshType(pos, normal, uv, color, indices);
              scene.addMesh(mesh);
            } // if (obj.name() == "Geometry")
          }
        } // if (section.name() == "Objects")
      }
      return true;
    }

  private:

    static size_t map(fbx_decoder::Mapping m, size_t polygon_index, size_t pvi, size_t vi) {
      switch (m) {
        case Mapping::ByPolygon: return polygon_index;
        case Mapping::ByPolygonVertex: return pvi;
        case Mapping::ByVertex: return vi;
        case Mapping::ByEdge: return 0; // todo
        default:
        case Mapping::AllSame: return 0;
      }
    }

    void init(const char *begin, const char *end) {
      begin_ = begin;
      end_ = end;

      if (end < begin + 27+4 || memcmp(begin, "Kaydara FBX Binary  ", 20)) {
        bad_fbx();
      }

      const char *p = begin + 23;
      std::string text;

      std::vector <const char *> ends;
      std::uint32_t version = u4(p);
      p += 4;

      while (u4(p)) {
        p = begin + u4(p);
      }
      end_offset = p - begin;
    }

    void dump(std::ostream &os, node &n, int depth, char *tmp, size_t tmp_size) const {
      snprintf(tmp, tmp_size, "%*sbegin(\"%s\");\n", depth*2, "", n.name().c_str());
      os << tmp;
      std::vector<uint64_t> ldata;
      std::vector<uint32_t> idata;
      for (auto p : n.get_props()) {
        snprintf(tmp, tmp_size, "%*s  %c(%s);\n", depth*2, "", p.kind(), ((std::string)p).c_str());
        os << tmp;
      }
      for (auto child : n) {
        dump(os, child, depth+1, tmp, tmp_size);
      }
      snprintf(tmp, tmp_size, "%*send(\"%s\");\n", depth*2, "", n.name().c_str());
      os << tmp;
    }

    friend std::ostream &operator<<(std::ostream &os, const fbx_decoder &fbx);

    void bad_fbx() { throw std::runtime_error("bad fbx"); }

    //std::vector<std::uint8_t> bytes;
    //scene the_scene;

    size_t end_offset;
    const char *begin_;
    const char *end_;
  };

  inline std::ostream &operator<<(std::ostream &os, const fbx_decoder &fbx) {
    std::vector<char> tmp(65536);
    for (auto p : fbx) {
      fbx.dump(os, p, 0, tmp.data(), tmp.size());
    }
    return os;
  }

}

#endif
