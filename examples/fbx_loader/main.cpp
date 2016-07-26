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

    std::ofstream txt1("1.txt", std::ios_base::binary);
    std::ofstream txt2("2.txt", std::ios_base::binary);

    const char *out_filename = "test.fbx";

    txt1 << fbx;

    //std::ofstream of(out_filename, std::ios_base::binary);
    meshutils::fbx_encoder encoder;
    //of.write((char*)encoder.bytes().data(), encoder.bytes().size());

    const char *beg = (char*)encoder.bytes().data();
    const char *end = beg + encoder.bytes().size();

    for (size_t i = 27; i != 1000; ++i) {
      printf("[%02x %02x%s", text[i] & 0xff, beg[i] & 0xff, text[i] != beg[i] ? "]\n" : "]");
    }
    printf("\n");

    meshutils::fbx_decoder fbx2(beg, end);
    txt2 << fbx2;
    
    return 0;
  } else {
    std::cerr << "uable to open file " << filename << "\n";
    return 1;
  }
  
}

