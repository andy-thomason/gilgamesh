
#include <meshutils/mesh.hpp>
#include <meshutils/box.hpp>
#include <meshutils/encoders/fbx_encoder.hpp>
#include <fstream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

int main() {
  meshutils::simple_mesh mesh;

  printf("building a cube\n");
  meshutils::box box{};
  box.build(mesh);

  printf("writing a fbx file\n");
  meshutils::fbx_encoder().saveMesh(mesh, "cube.fbx");
}
