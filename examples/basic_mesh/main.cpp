
#include <meshutils/mesh.hpp>
#include <meshutils/shapes/box.hpp>
#include <meshutils/shapes/sphere.hpp>
#include <meshutils/shapes/cylinder.hpp>
#include <meshutils/encoders/fbx_encoder.hpp>
#include <meshutils/decoders/fbx_decoder.hpp>
#include <fstream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

int main() {
  meshutils::color_mesh mesh;
  meshutils::fbx_encoder encoder;
  meshutils::fbx_decoder decoder;
  glm::mat4 mat;

  printf("building a cube\n");
  meshutils::box box{};
  box.build(mesh, mat, glm::vec4(1, 0, 0, 1));

  printf("writing cube.fbx\n");
  encoder.saveMesh(mesh, "cube.fbx");

  mesh.clear();

  printf("building a sphere\n");
  meshutils::sphere sphere{};
  sphere.build(mesh, mat, glm::vec4(0, 1, 0, 1));

  printf("writing sphere.fbx\n");
  encoder.saveMesh(mesh, "sphere.fbx");

  mesh.clear();

  printf("building a cylinder\n");
  meshutils::cylinder cylinder{};
  cylinder.build(mesh, mat, glm::vec4(0, 0, 1, 1));

  printf("writing cylinder.fbx\n");
  encoder.saveMesh(mesh, "cylinder.fbx");

  mesh.clear();

  printf("building a composite object\n");

  box.build(mesh, mat, glm::vec4(1, 0, 0, 1));
  mat[3].x = -2;
  sphere.build(mesh, mat, glm::vec4(0, 1, 0, 1));
  mat[3].x = 2;
  cylinder.build(mesh, mat, glm::vec4(0, 0, 1, 1));

  printf("writing composite.fbx\n");
  encoder.saveMesh(mesh, "composite.fbx");
}
