////////////////////////////////////////////////////////////////////////////////
//
// (C) Andy Thomason 2016
//
// example of FBX file decoder and encoder
//
////////////////////////////////////////////////////////////////////////////////

#include <gilgamesh/mesh.hpp>
#include <gilgamesh/scene.hpp>
#include <gilgamesh/decoders/fbx_decoder.hpp>
#include <gilgamesh/encoders/fbx_encoder.hpp>

int main(int argc, char **argv) {
  const char *filename = nullptr;
  for (int i = 1; i != argc; ++i) {
    const char *arg = argv[i];
    if (arg[0] == '-') {
      printf("invalid argument %s\n", arg);
    } else {
      filename = arg;
    }
  }

  if (filename == nullptr) {
    filename = CMAKE_SOURCE "/examples/data/cube.fbx";
  }

  gilgamesh::fbx_decoder decoder;

  gilgamesh::scene scene;
  if (!decoder.loadScene<gilgamesh::color_mesh>(scene, filename)) {
    std::cerr << "uable to open file " << filename << "\n";
    return 1;
  }

  gilgamesh::fbx_encoder encoder;
  auto bytes = encoder.saveScene(scene);

  std::ofstream("out.fbx", std::ios::binary).write((char*)bytes.data(), bytes.size());
}

