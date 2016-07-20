
#include <meshutils/mesh.hpp>
#include <fstream>

int main() {
  meshutils::simple_mesh mesh;

  mesh.addCube();

  std::ofstream of("mesh.csv");
  mesh.writeCSV(of);
}
