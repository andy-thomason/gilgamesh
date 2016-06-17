
#include <meshutils/mesh.hpp>

int main() {
  meshutils::simple_mesh mesh;

  mesh.addCube();

  std::ofstream of("mesh.csv");
  mesh.writeCSV(of);
}
