////////////////////////////////////////////////////////////////////////////////
//
// (C) Andy Thomason 2016
//
// High performance Protien Data Bank file format reader
// 
// PDB files are Fortran-style text files containing positions of atoms in molecules.

#ifndef MESHUTILS_pdb_decoder_INCLUDED
#define MESHUTILS_pdb_decoder_INCLUDED

#include <iostream>
#include <cstdint>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>

#include <glm/glm.hpp>


// https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)
// https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/framepdbintro.html
namespace meshutils {
  class pdb_decoder {
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
      std::string atomName() const { return std::string(p_ - 1 + 13, p_ + 16); }
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

    pdb_decoder(const uint8_t *begin, const uint8_t *end) {
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
                hetAtoms_.emplace_back(p, eol);
              }
            } break;
            case 'C': {
              if (p + 5 < eol && !memcmp(p, "CONECT", 6)) {
                // COLUMNS       DATA  TYPE      FIELD        DEFINITION
                // -------------------------------------------------------------------------
                //  1 -  6        Record name    "CONECT"
                //  7 - 11       Integer        serial       Atom  serial number
                //  12 - 16        Integer        serial       Serial number of bonded atom
                //  17 - 21        Integer        serial       Serial  number of bonded atom
                //  22 - 26        Integer        serial       Serial number of bonded atom
                //  27 - 31        Integer        serial       Serial number of bonded atom
                int a0 = atoi(p + 1 + 7, p + 11);
                int a1 = atoi(p + 1 + 12, p + 16);
                int a2 = atoi(p + 1 + 17, p + 21);
                int a3 = atoi(p + 1 + 22, p + 26);
                int a4 = atoi(p + 1 + 27, p + 31);
                if (a0 && a1) connections_.emplace_back(a0, a1);
                if (a0 && a2) connections_.emplace_back(a0, a2);
                if (a0 && a3) connections_.emplace_back(a0, a3);
                if (a0 && a4) connections_.emplace_back(a0, a4);
              }
            } break;
          }
        }
        p = next_p;
      }
    }

    const std::vector<atom> &atoms() const { return atoms_; }
    const std::vector<atom> &hetAtoms() const { return hetAtoms_; }

    // return the set of chains used in this PDB file (ie. "ABCD")
    std::string chains() const {
      bool used[256] = {};
      for (auto &p : atoms_) {
        used[p.chainID()] = true;
      }
      std::string result;
      for (size_t i = 32; i != 127; ++i) {
        if (used[i]) result.push_back((char)i);
      }
      return std::move(result);
    }

    // return a vector of positions filtered by chain id
    std::vector<glm::vec3> pos(char chainID = '?') const {
      std::vector<glm::vec3> result;
      for (auto &p : atoms_) {
        if (chainID == '?' || p.chainID() == chainID) {
          result.push_back(glm::vec3(p.x(), p.y(), p.z()));
        }
      }
      return std::move(result);
    }

    // Van Der Walls radii of atoms
    std::vector<float> radii(char chainID = '?') const {
      std::vector<float> result;
      for (auto &p : atoms_) {
        if (chainID == '?' || p.chainID() == chainID) {
          result.push_back(vdvRadius(p.element()));
        }
      }
      return std::move(result);
    }

    // Return colours for terminal atoms of amino acids.
    std::vector<glm::vec4> colorsByFunction(char chainID = '?') const {
      std::vector<glm::vec4> result;
      for (auto &p : atoms_) {
        if (chainID == '?' || p.chainID() == chainID) {
          std::string atom = p.atomName();
          std::string resName = p.resName();
          if (
            (atom == " NZ " && resName == "LYS") ||
            (atom == " NH1" && resName == "ARG") ||
            (atom == " NH2" && resName == "ARG") ||
            (atom == " ND1" && resName == "HIS") ||
            (atom == " NE2" && resName == "HIS")
          ) {
            // Positive: blue
            result.push_back(glm::vec4(0, 0, 1, 1));
          } else if (
            (atom == " OE1" && resName == "GLU") ||
            (atom == " OE2" && resName == "GLU") ||
            (atom == " OD1" && resName == "ASP") ||
            (atom == " OD2" && resName == "ASP")
          ) {
            // Negative: red
            result.push_back(glm::vec4(1, 0, 0, 1));
          } else {
            // default: white
            result.push_back(glm::vec4(1, 1, 1, 1));
          }
        }
      }
      return std::move(result);
    }

    // return implicit connections between atoms
    std::vector<std::pair<int, int> > implicitConnections(char chainID = '?') const {
      std::vector<std::pair<int, int> > result;

      // see http://www.bmrb.wisc.edu/referenc/commonaa.php?asp
      int prevC = -1;
      int idx = 0;
      for (auto &p : atoms_) {
        if (chainID == '?' || p.chainID() == chainID) {
          std::string atom = p.atomName();
          if (atom == " N  ") {
            // backbone
            int N = idx, CA = idx+1, C = idx+2, O = idx+3, CB = idx + 4;
            if (prevC != -1) result.emplace_back(prevC, N);
            result.emplace_back(N, CA);
            result.emplace_back(CA, C);
            result.emplace_back(CA, O);

            std::string resName = p.resName();
            // CB
            if (resName != "GLY") {
              result.emplace_back(CA, CB);
            }

            //I can't guarantee these are correct right now.
            if (resName == "ASP") {
              int CG = CB + 1;
              result.emplace_back(CB, CG);
              int OD1 = CG + 1;
              int OD2 = CG + 2;
              result.emplace_back(CG, OD1);
              result.emplace_back(CG, OD2);
            //} else if (resName == "ALA") {
            } else if (resName == "CYS") {
              int SG = idx + 4;
              result.emplace_back(CB, SG);
            } else if (resName == "GLU") {
              int CG = CB + 1;
              result.emplace_back(CB, CG);
              int CD = CG + 1;
              result.emplace_back(CG, CD);
              int OE1 = CD + 1;
              int OE2 = CD + 2;
              result.emplace_back(CD, OE1);
              result.emplace_back(CD, OE2);
            } else if (resName == "PHE") {
              int CG = CB + 1;
              result.emplace_back(CB, CG);
              int CD1 = CG + 1;
              int CD2 = CG + 2;
              result.emplace_back(CG, CD1);
              result.emplace_back(CG, CD2);
              int CE1 = CG + 3;
              int CE2 = CG + 4;
              result.emplace_back(CD1, CE1);
              result.emplace_back(CD2, CE2);
              int CZ = CG + 5;
              result.emplace_back(CE1, CZ);
              result.emplace_back(CE2, CZ);
            //} else if (resName == "GLY") {
            } else if (resName == "HIS") {
              int CG = CB + 1;
              result.emplace_back(CB, CG);
              int ND1 = CG + 1;
              int CD2 = CG + 2;
              result.emplace_back(CG, ND1);
              result.emplace_back(CG, CD2);
              int CE1 = CG + 3;
              int NE2 = CG + 4;
              result.emplace_back(ND1, CE1);
              result.emplace_back(CD2, NE2);
              result.emplace_back(CE1, NE2);
            } else if (resName == "ILE") {
              int CG1 = CB + 1;
              int CG2 = CB + 2;
              result.emplace_back(CB, CG1);
              result.emplace_back(CB, CG2);
              result.emplace_back(CG1, CG2);
            } else if (resName == "LYS") {
              int CG = CB + 1;
              result.emplace_back(CB, CG);
              int CD = CG + 1;
              result.emplace_back(CG, CD);
              int CE = CD + 1;
              result.emplace_back(CD, CE);
              int NZ = CE + 1;
              result.emplace_back(CE, NZ);
            } else if (resName == "LEU") {
              int CG = CB + 1;
              result.emplace_back(CB, CG);
              int CD1 = CG + 1;
              int CD2 = CG + 2;
              result.emplace_back(CG, CD1);
              result.emplace_back(CG, CD2);
            } else if (resName == "MET") {
              int CG = CB + 1;
              result.emplace_back(CB, CG);
              int SD = CG + 1;
              result.emplace_back(CG, SD);
              int CE = SD + 1;
              result.emplace_back(SD, CE);
            } else if (resName == "ASN") {
              int CG = CB + 1;
              result.emplace_back(CB, CG);
              int OD1 = CG + 1;
              int ND2 = CG + 2;
              result.emplace_back(CG, OD1);
              result.emplace_back(CG, ND2);
            } else if (resName == "PRO") {
              int CG = CB + 1;
              result.emplace_back(CB, CG);
              int CD = CG + 1;
              result.emplace_back(CG, CD);
            } else if (resName == "GLN") {
              int CG = CB + 1;
              result.emplace_back(CB, CG);
              int CD = CG + 1;
              result.emplace_back(CG, CD);
              int OE1 = CD + 1;
              int NE2 = CD + 2;
              result.emplace_back(CD, OE1);
              result.emplace_back(CD, NE2);
            } else if (resName == "ARG") {
              int CG = CB + 1;
              result.emplace_back(CB, CG);
              int CD = CG + 1;
              result.emplace_back(CG, CD);
              int NE = CD + 1;
              result.emplace_back(CD, NE);
              int CZ = NE + 1;
              result.emplace_back(NE, CZ);
              int NH1 = CZ + 1;
              int NH2 = CZ + 2;
              result.emplace_back(CZ, NH1);
              result.emplace_back(CZ, NH2);
            } else if (resName == "SER") {
              int OG = CB + 1;
              result.emplace_back(CB, OG);
            } else if (resName == "THR") {
              int OG1 = CB + 1;
              int CG2 = CB + 2;
              result.emplace_back(CB, OG1);
              result.emplace_back(CB, CG2);
            } else if (resName == "VAL") {
              int CG1 = CB + 1;
              int CG2 = CB + 2;
              result.emplace_back(CB, CG1);
              result.emplace_back(CB, CG2);
            } else if (resName == "TRP") {
              int CG = CB + 1;
              result.emplace_back(CB, CG);
              int CD1 = CG + 1;
              int CD2 = CG + 2;
              result.emplace_back(CG, CD1);
              result.emplace_back(CG, CD2);
              int NE1 = CG + 3;
              int CE2 = CG + 4;
              int CE3 = CG + 5;
              result.emplace_back(CD1, NE1);
              result.emplace_back(CD2, CE3);
              result.emplace_back(NE1, CE2);
              int CZ2 = CG + 5;
              int CZ3 = CG + 6;
              result.emplace_back(CE2, CZ2);
              result.emplace_back(CE3, CZ3);
              int CH2 = CG + 7;
              result.emplace_back(CZ2, CH2);
              result.emplace_back(CZ3, CH2);
            } else if (resName == "TYR") {
              int CG = CB + 1;
              result.emplace_back(CB, CG);
              int CD1 = CG + 1;
              int CD2 = CG + 2;
              result.emplace_back(CG, CD1);
              result.emplace_back(CG, CD2);
              int CE1 = CG + 3;
              int CE2 = CG + 4;
              result.emplace_back(CD1, CE1);
              result.emplace_back(CD2, CE2);
              int CZ = CG + 5;
              result.emplace_back(CE1, CZ);
              result.emplace_back(CE2, CZ);
              int OH = CZ + 1;
              result.emplace_back(CZ, OH);
            }

          }
        }
      }
      return result;
    }
  private:
    std::vector<atom> atoms_;
    std::vector<atom> hetAtoms_;
    std::vector<std::pair<int, int> > connections_;

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
