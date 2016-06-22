////////////////////////////////////////////////////////////////////////////////
//
// (C) Andy Thomason 2016
//
// High performance Protien Data Bank file format reader
// 
// PDB files are Fortran-style text files containing positions of atoms in molecules.

#ifndef MESHUTILS_PDB_FILE_INCLUDED
#define MESHUTILS_PDB_FILE_INCLUDED

#include <iostream>
#include <cstdint>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <glm/glm.hpp>

// https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)
namespace meshutils {
  class pdb_file {
  public:
    class atom {
      const uint8_t *p_;
      const uint8_t *eol_;
    public:
      atom(const uint8_t *p, const uint8_t *eol) : p_(p), eol_(eol) {
      }

      // http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
      //  7 - 11        Integer       serial       Atom  serial number.
      // 13 - 16        Atom          name         Atom name.
      // 17             Character     altLoc       Alternate location indicator.
      // 18 - 20        Residue name  resName      Residue name.
      // 22             Character     chainID      Chain identifier.
      // 23 - 26        Integer       resSeq       Residue sequence number.
      // 27             AChar         iCode        Code for insertion of residues.
      // 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
      // 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
      // 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
      // 55 - 60        Real(6.2)     occupancy    Occupancy.
      // 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
      // 77 - 78        LString(2)    element      Element symbol, right-justified.
      // 79 - 80        LString(2)    charge       Charge  on the atom.
      int serial() const { return atoi(p_ - 1 + 7, p_ + 11); }
      std::string name() const { return std::string(p_ - 1 + 13, p_ + 16); }
      char altLoc() const { return (char)p_[-1+17]; }
      std::string resName() const { return std::string(p_ - 1 + 18, p_ + 20); }
      char chainID() const { return (char)p_[-1+22]; }
      int resSeq() const { return atoi(p_ - 1 + 23, p_ + 26); }
      char iCode() const { return (char)p_[-1+27]; }
      float x() const { return atof(p_ - 1 + 31, p_ + 38); }
      float y() const { return atof(p_ - 1 + 39, p_ + 46); }
      float z() const { return atof(p_ - 1 + 47, p_ + 54); }
      float occupancy() const { return atof(p_ - 1 + 55, p_ + 60); }
      float tempFactor() const { return atof(p_ - 1 + 61, p_ + 66); }
      std::string element() const { return std::string(p_ - 1 + 77, p_ + 78); }
      std::string charge() const { return std::string(p_ - 1 + 79, p_ + 80); }
    };

    pdb_file(const uint8_t *begin, const uint8_t *end) {
      for (const uint8_t *p = begin; p != end; ) {
        const uint8_t *eol = p;
        while (eol != end && *eol != '\n') ++eol;
        const uint8_t *next_p = eol != end ? eol + 1 : end;
        while (eol != p && (*eol == '\r' || *eol == '\n')) --eol;
        if (p != eol) {
          switch (*p) {
            case 'A': {
              if (p + 5 < eol && !memcmp(p, "ATOM  ", 6)) {
                atoms_.emplace_back(p, eol);
              }
            } break;
            case 'H': {
              if (p + 5 < eol && !memcmp(p, "HETATM", 6)) {
                hetatoms_.emplace_back(p, eol);
              }
            } break;
          }
        }
        p = next_p;
      }
    }

    const std::vector<atom> &atoms() const { return atoms_; }

    std::string chains() const {
      bool used[256] = {};
      for (auto &p : atoms_) {
        used[p.chainID()] = true;
      }
      std::string result;
      for (size_t i = 0; i != 128; ++i) {
        if (used[i]) result.push_back((char)i);
      }
      return result;
    }

    std::vector<glm::vec3> pos(char chainID = '?') const {
      std::vector<glm::vec3> result;
      for (auto &p : atoms_) {
        if (chainID == '?' || p.chainID() == chainID) {
          result.push_back(glm::vec3(p.x(), p.y(), p.z()));
        }
      }
      return result;
    }

    // Van Der Walls radii of atoms
    std::vector<float> radii(char chainID = '?') const {
      std::vector<float> result;
      for (auto &p : atoms_) {
        if (chainID == '?' || p.chainID() == chainID) {
          result.push_back(vdvRadius(p.name()));
        }
      }
      return result;
    }
  public:
  private:
    std::vector<atom> atoms_;
    std::vector<atom> hetatoms_;

    static int atoi(const uint8_t *b, const uint8_t *e) {
      while (b != e && *b == ' ') ++b;
      int n = 0;
      while (b != e && *b >= '0' && *b <= '9') n = n * 10 + *b++ - '0';
      return n;
    }

    static float atof(const uint8_t *b, const uint8_t *e) {
      while (b != e && *b == ' ') ++b;
      float n = 0;
      float s = 1.0f;
      if (b != e && *b == '-') { s = -s; b++; }
      while (b != e && *b >= '0' && *b <= '9') n = n * 10 + *b++ - '0';
      if (b != e && *b == '.') {
        ++b;
        float frac = 0, p10 = 1;
        while (b != e && *b >= '0' && *b <= '9') p10 *= 10, frac = frac * 10 + *b++ - '0';
        n += frac / p10;
      }
      if (b != e && (*b == 'e'||*b == 'E')) {
        ++b;
        int es = 1;
        if (b != e && *b == '-') { es = -es; b++; }
        int exp = 0;
        while (b != e && *b >= '0' && *b <= '9') exp = exp * 10 + *b++ - '0';
        return s * n * std::pow(10.0f, exp * es);
      } else {
        return s * n;
      }
    }

    static float vdvRadius(const std::string &element) {
      // https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)

      struct data_t { char name[4]; short vdv; };

      static const data_t data[] = {
        "H ", 120, "C ", 170, "N ", 155, "O ", 152, "S ", 180, "HE", 140, "LI", 182, "BE", 153, "B ", 192, 
        "F ", 147, "NE", 154, "NA", 227, "MG", 173, "AL", 184, "SI", 210, "P ", 180,
        "CL", 175, "AR", 188, "K ", 275, "CA", 231, "SC", 211, "NI", 163, "CU", 140, "ZN", 139,
        "GA", 187, "GE", 211, "AS", 185, "SE", 190, "BR", 185, "KR", 202, "RB", 303, "SR", 249, 
        "PD", 163, "AG", 172, "CD", 158, "IN", 193, "SN", 217, "SB", 206, "TE", 206, "I ", 198, 
        "XE", 216, "CS", 343, "BA", 268, "PT", 175, "AU", 166, "HG", 155, "TL", 196, "PB", 202,
        "BI", 207, "PO", 197, "AT", 202, "RN", 220, "FR", 348, "RA", 283, "U ", 186,
      };

      char e0 = element[1];
      char e1 = element[2];
      for (const data_t &d : data) {
        if (e0 == d.name[0] && e1 == d.name[1]) {
          return d.vdv * 0.01f;
        }
      }
      return 1.2f;
    }
  };
}

#endif
