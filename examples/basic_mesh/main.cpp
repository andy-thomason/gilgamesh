
#include <meshutils/mesh.hpp>
#include <meshutils/encoders/ply_encoder.hpp>
#include <fstream>

int main() {
  meshutils::simple_mesh mesh;

  printf("building a cube\n");
  mesh.addCube();

  printf("writing a stanford ply file\n");
  meshutils::ply_encoder enc;
  std::ofstream of("cube.ply", std::ios::binary);
  enc.encode(mesh, of);
}
