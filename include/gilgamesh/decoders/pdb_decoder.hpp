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
namespace gilgamesh {
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

      bool resNameIs(const char *name) const { return p_[17] == name[0] && p_[18] == name[1] && p_[19] == name[2]; }
      bool atomNameIs(const char *name) const { return p_[12] == name[0] && p_[13] == name[1] && p_[14] == name[2] && p_[15] == name[3]; }
      bool elementIs(const char *name) const { return p_[76] == name[0] && p_[77] == name[1]; }
      bool chargeIs(const char *name) const { return p_[78] == name[0] && p_[79] == name[1]; }

      glm::vec4 colorByFunction() const {
        std::string atom = atomName();
        if (
          (atomNameIs(" NZ ") && resNameIs("LYS")) ||
          (atomNameIs(" NH1") && resNameIs("ARG")) ||
          (atomNameIs(" NH2") && resNameIs("ARG")) ||
          (atomNameIs(" ND1") && resNameIs("HIS")) ||
          (atomNameIs(" NE2") && resNameIs("HIS"))
        ) {
          // Positive: blue
          return glm::vec4(0, 0, 1, 1);
        } else if (
          (atomNameIs(" OE1") && resNameIs("GLU")) ||
          (atomNameIs(" OE2") && resNameIs("GLU")) ||
          (atomNameIs(" OD1") && resNameIs("ASP")) ||
          (atomNameIs(" OD2") && resNameIs("ASP"))
        ) {
          // Negative: red
          return glm::vec4(1, 0, 0, 1);
        } else {
          return glm::vec4(1, 1, 1, 1);
        }
      }

      // for the start of an amino acid, add implicit connections 
      void addImplicitConnections(std::vector<std::pair<int, int> > &result, int idx) const {

        // backbone
        int N = idx, CA = idx+1, C = idx+2, O = idx+3, CB = idx + 4;
        result.emplace_back(N, CA);
        result.emplace_back(CA, C);
        result.emplace_back(CA, O);

        if (!resNameIs("GLY")) {
          result.emplace_back(CA, CB);
        }

        // I can't guarantee these are correct right now.
        if (resNameIs("ASP")) {
          int CG = CB + 1;
          result.emplace_back(CB, CG);
          int OD1 = CG + 1;
          int OD2 = CG + 2;
          result.emplace_back(CG, OD1);
          result.emplace_back(CG, OD2);
        } else if (resNameIs("ALA")) {
        } else if (resNameIs("CYS")) {
          int SG = idx + 4;
          result.emplace_back(CB, SG);
        } else if (resNameIs("GLU")) {
          int CG = CB + 1;
          result.emplace_back(CB, CG);
          int CD = CG + 1;
          result.emplace_back(CG, CD);
          int OE1 = CD + 1;
          int OE2 = CD + 2;
          result.emplace_back(CD, OE1);
          result.emplace_back(CD, OE2);
        } else if (resNameIs("PHE")) {
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
        } else if (resNameIs("GLY")) {
        } else if (resNameIs("HIS")) {
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
        } else if (resNameIs("ILE")) {
          int CG1 = CB + 1;
          int CG2 = CB + 2;
          result.emplace_back(CB, CG1);
          result.emplace_back(CB, CG2);
          result.emplace_back(CG1, CG2);
        } else if (resNameIs("LYS")) {
          int CG = CB + 1;
          result.emplace_back(CB, CG);
          int CD = CG + 1;
          result.emplace_back(CG, CD);
          int CE = CD + 1;
          result.emplace_back(CD, CE);
          int NZ = CE + 1;
          result.emplace_back(CE, NZ);
        } else if (resNameIs("LEU")) {
          int CG = CB + 1;
          result.emplace_back(CB, CG);
          int CD1 = CG + 1;
          int CD2 = CG + 2;
          result.emplace_back(CG, CD1);
          result.emplace_back(CG, CD2);
        } else if (resNameIs("MET")) {
          int CG = CB + 1;
          result.emplace_back(CB, CG);
          int SD = CG + 1;
          result.emplace_back(CG, SD);
          int CE = SD + 1;
          result.emplace_back(SD, CE);
        } else if (resNameIs("ASN")) {
          int CG = CB + 1;
          result.emplace_back(CB, CG);
          int OD1 = CG + 1;
          int ND2 = CG + 2;
          result.emplace_back(CG, OD1);
          result.emplace_back(CG, ND2);
        } else if (resNameIs("PRO")) {
          int CG = CB + 1;
          result.emplace_back(CB, CG);
          int CD = CG + 1;
          result.emplace_back(CG, CD);
        } else if (resNameIs("GLN")) {
          int CG = CB + 1;
          result.emplace_back(CB, CG);
          int CD = CG + 1;
          result.emplace_back(CG, CD);
          int OE1 = CD + 1;
          int NE2 = CD + 2;
          result.emplace_back(CD, OE1);
          result.emplace_back(CD, NE2);
        } else if (resNameIs("ARG")) {
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
        } else if (resNameIs("SER")) {
          int OG = CB + 1;
          result.emplace_back(CB, OG);
        } else if (resNameIs("THR")) {
          int OG1 = CB + 1;
          int CG2 = CB + 2;
          result.emplace_back(CB, OG1);
          result.emplace_back(CB, CG2);
        } else if (resNameIs("VAL")) {
          int CG1 = CB + 1;
          int CG2 = CB + 2;
          result.emplace_back(CB, CG1);
          result.emplace_back(CB, CG2);
        } else if (resNameIs("TRP")) {
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
        } else if (resNameIs("TYR")) {
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

      glm::vec4 colorByElement() const {
        // https://en.wikipedia.org/wiki/CPK_coloring

        struct data_t { char name[4]; uint32_t color; };

        static const data_t cpk[] = {
          " H", 0xffffff, " C", 0x222222, " N", 0x2233ff, " O", 0xff2200, " S", 0xdddd00,
        };
        static const data_t jmol[] = {
          " H", 0xffffff, " C", 0x909090, " N", 0x3050F8, " O", 0xFF0D0D, " F", 0x90E050,
          "NA", 0xAB5CF2, "MG", 0x8AFF00, "AL", 0xBFA6A6, "SI", 0xF0C8A0,
          " P", 0xFF8000, " S", 0xFFFF30, "CL", 0x1FF01F, "AR", 0x80D1E3,
          " K", 0x8F40D4, "CA", 0x3DFF00,
        };

        char e0 = p_[76];
        char e1 = p_[77];
        uint32_t color = 0xdd77ff;
        for (const data_t &d : jmol) {
          if (e0 == d.name[0] && e1 == d.name[1]) {
            color = d.color;
            break;
          }
        }
        return glm::vec4((color >> 16) * (1.0f/255), ((color >> 8)&0xff) * (1.0f/255), (color & 0xff) * (1.0f/255), 1.0f);
      }

      float vanDerVaalsRadius() const {
        // https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)

        struct data_t { char name[4]; short vdv; };

        static const data_t data[] = {
          " H", 120, " C", 170, " N", 155, " O", 152, " S", 180, "HE", 140, "LI", 182, "BE", 153, " B", 192, 
          " F", 147, "NE", 154, "NA", 227, "MG", 173, "AL", 184, "SI", 210, " P", 180,
          "CL", 175, "AR", 188, " K", 275, "CA", 231, "SC", 211, "NI", 163, "CU", 140, "ZN", 139,
          "GA", 187, "GE", 211, "AS", 185, "SE", 190, "BR", 185, "KR", 202, "RB", 303, "SR", 249, 
          "PD", 163, "AG", 172, "CD", 158, "IN", 193, "SN", 217, "SB", 206, "TE", 206, " I", 198, 
          "XE", 216, "CS", 343, "BA", 268, "PT", 175, "AU", 166, "HG", 155, "TL", 196, "PB", 202,
          "BI", 207, "PO", 197, "AT", 202, "RN", 220, "FR", 348, "RA", 283, " U", 186,
        };

        char e0 = p_[76];
        char e1 = p_[77];
        for (const data_t &d : data) {
          if (e0 == d.name[0] && e1 == d.name[1]) {
            return d.vdv * 0.01f;
          }
        }
        return 1.2f;
      }
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

    static int addImplicitConnections(std::vector<std::pair<int, int> > out, const atom *atomb, const atom *atome, int N_idx) {
      static const char table[][5] = {
        "ASP",
          " CB ", " CG ",
          " CG ", " OD1",
          " CG ", " OD2",
        "ALA",
        "CYS",
          " CB ", " SG ",
        "GLU",
          " CB ", " CG ",
          " CG ", " CD ",
          " CD ", " OE1",
          " CD ", " OE2",
        "PHE",
          " CB ", " CG ",
          " CG ", " CD1",
          " CG ", " CD2",
          " CD1", " CE1",
          " CD2", " CE2",
          " CE1", " CZ ",
          " CE2", " CZ ",
        "GLY",
        "HIS",
          " CB ", " CG ",
          " CG ", " ND1",
          " CG ", " CD2",
          " ND1", " CE1",
          " CD2", " NE2",
          " CE1", " NE2",
        "ILE",
          " CB ", " CG1",
          " CB ", " CG2",
          " CG1", " CG2",
        "LYS",
          " CB ", " CG ",
          " CG ", " CD ",
          " CD ", " CE ",
          " CE ", " NZ ",
        "LEU",
          " CB ", " CG ",
          " CG ", " CD1",
          " CG ", " CD2",
        "MET",
          " CB ", " CG ",
          " CG ", " SD ",
          " SD ", " CE ",
        "ASN",
          " CB ", " CG ",
          " CG ", " OD1",
          " CG ", " ND2",
        "PRO",
          " CB ", " CG ",
          " CG ", " CD ",
        "GLN",
          " CB ", " CG ",
          " CG ", " CD ",
          " CD ", " OE1",
          " CD ", " NE2",
        "ARG",
          " CB ", " CG ",
          " CG ", " CD ",
          " CD ", " NE ",
          " NE ", " CZ ",
          " CZ ", " NH1",
          " CZ ", " NH2",
        "SER",
          " CB ", " OG ",
        "THR",
          " CB ", " OG1",
          " CB ", " CG2",
        "VAL",
          " CB ", " CG1",
          " CB ", " CG2",
        "TRP",
          " CB ", " CG ",
          " CG ", " CD1",
          " CG ", " CD2",
          " CD1", " NE1",
          " CD2", " CE3",
          " NE1", " CE2",
          " CE2", " CZ2",
          " CE3", " CZ3",
          " CZ2", " CH2",
          " CZ3", " CH2",
        "TYR",
          " CB ", " CG ",
          " CG ", " CD1",
          " CG ", " CD2",
          " CD1", " CE1",
          " CD2", " CE2",
          " CE1", " CZ ",
          " CE2", " CZ ",
          " CZ ", " OH ",
        ""
      };

      // find an atom relative to "N" in an amino acid.
      auto find_atom = [atomb, atome, N_idx](const char *name) {
        for (auto p = atomb; p != atome; ++p) {
          if (p->atomNameIs(name)) {
            return int(p - atomb) + N_idx;
          }
        }
        return -1;
      };


      int C_idx = find_atom(" C  ");
      int O_idx = find_atom(" O  ");
      int CA_idx = find_atom(" CA ");
      int CB_idx = find_atom(" CB ");

      //printf("find %s N%d C%d O%d CA%d CB%d\n", atomb->resName().c_str(), N_idx, C_idx, O_idx, CA_idx, CB_idx);

      out.emplace_back(N_idx, C_idx);
      if (CA_idx != -1) out.emplace_back(C_idx, CA_idx);
      if (CA_idx != -1 && CB_idx != -1) out.emplace_back(CA_idx, CB_idx);
      out.emplace_back(C_idx, O_idx);

      for (size_t i = 0; table[i][0]; ++i) {
        if (table[i][0] >= 'A' && atomb->resNameIs(table[i])) {
          //printf("%s\n", table[i]);
          ++i;
          while (table[i][0] == ' ') {
            //printf("  %s %s\n", table[i], table[i+1]);
            int from = find_atom(table[i]);
            int to = find_atom(table[i+1]);
            i += 2;
            //printf("  %d..%d\n", from, to);
            out.emplace_back(from, to);
            //if (from == -1 || to == -1) printf("OUCH!\n");
          }
          break;
        }
      }

      return C_idx;
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

  };
}

#endif
