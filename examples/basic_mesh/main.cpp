
#include <gilgamesh/mesh.hpp>
#include <gilgamesh/shapes/box.hpp>
#include <gilgamesh/shapes/sphere.hpp>
#include <gilgamesh/shapes/cylinder.hpp>
#include <gilgamesh/encoders/fbx_encoder.hpp>
#include <gilgamesh/decoders/fbx_decoder.hpp>
#include <fstream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

int main() {
  gilgamesh::color_mesh mesh;
  gilgamesh::fbx_encoder encoder;
  gilgamesh::fbx_decoder decoder;
  glm::mat4 mat;

  printf("building a cube\n");
  gilgamesh::box box{};
  box.build(mesh, mat, glm::vec4(1, 0, 0, 1));

  printf("writing cube.fbx\n");
  encoder.saveMesh(mesh, "cube.fbx");

  mesh.clear();

  printf("building a sphere\n");
  gilgamesh::sphere sphere{};
  sphere.build(mesh, mat, glm::vec4(0, 1, 0, 1));

  printf("writing sphere.fbx\n");
  encoder.saveMesh(mesh, "sphere.fbx");

  mesh.clear();

  printf("building a cylinder\n");
  gilgamesh::cylinder cylinder{};
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
