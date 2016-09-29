
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
  glm::mat4 mat;

  mesh.clear();
  gilgamesh::sphere sphere{};
  sphere.build(mesh, mat, glm::vec4(0, 1, 0, 1), 2, false);

  mesh.writeCSV("1.csv");
  sphere.build(mesh, mat, glm::vec4(0, 1, 0, 1), 2, true);

  mesh.writeCSV("2.csv");
}
