////////////////////////////////////////////////////////////////////////////////
//
// (C) Andy Thomason 2016
//
// example of FBX file decoder and encoder
//
////////////////////////////////////////////////////////////////////////////////

#include <meshutils/mesh.hpp>
#include <meshutils/scene.hpp>
#include <meshutils/decoders/fbx_decoder.hpp>
#include <meshutils/encoders/fbx_encoder.hpp>

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

  std::ifstream file(filename, std::ios_base::binary);
  std::vector<char> text;
  if (!file.eof() && !file.fail()) {
    file.seekg(0, std::ios_base::end);
    text.resize((size_t)file.tellg());

    file.seekg(0, std::ios_base::beg);
    file.read((char*)text.data(), text.size());
    meshutils::fbx_decoder fbx(text.data(), text.data() + text.size());

    // dump input to text file
    std::ofstream("1.txt", std::ios_base::binary) << fbx;

    meshutils::scene scene;
    fbx.loadScene<meshutils::color_mesh>(scene);
    std::cout << scene.meshes().size() << "\n";

    meshutils::fbx_encoder encoder;

    std::vector<uint8_t> bytes = encoder.saveScene(scene);

    const char *beg = (char*)bytes.data();
    const char *end = beg + bytes.size();

    // dump output to text file (so we can diff them)
    meshutils::fbx_decoder fbx2(beg, end);
    std::ofstream("2.txt", std::ios_base::binary) << fbx2;

    std::ofstream outfile("out.fbx", std::ios_base::binary);
    outfile.write(beg, end-beg);
    return 0;
  } else {
    std::cerr << "uable to open file " << filename << "\n";
    return 1;
  }
}

