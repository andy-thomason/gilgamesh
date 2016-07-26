////////////////////////////////////////////////////////////////////////////////
//
// (C) Andy Thomason 2016
//
////////////////////////////////////////////////////////////////////////////////

#include <meshutils/mesh.hpp>
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

    std::cout << fbx;

    const char *out_filename = "test.fbx";
    std::ofstream of(out_filename, std::ios_base::binary);
    meshutils::fbx_encoder encoder(of);
    return 0;
  } else {
    std::cerr << "uable to open file " << filename << "\n";
    return 1;
  }
  
}

