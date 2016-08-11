////////////////////////////////////////////////////////////////////////////////
//
// (C) Andy Thomason 2016
//
// example of FBX file decoder and encoder
//
////////////////////////////////////////////////////////////////////////////////

#include <meshutils/decoders/fbx_decoder.hpp>

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

  meshutils::fbx_decoder fbx;
  
  fbx.dump(std::cout, filename);
}

