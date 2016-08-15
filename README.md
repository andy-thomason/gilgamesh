Gilgamesh
=========

Mesh utilities for basic meshes, primitives, marching cubes, fbx files, lighting and CSG in modern C++

Gilgamesh is GPL-free so feel free to use it in your project, commercial or otherwise.

To get started, you will first need to clone this project and then:

```
git submodule init
git submodule update
```

This will import glm and minizip, required by this library.

Once you have the code, you will need CMake to build it.

Example:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
```

This builds either a Makefile (Linux/Mac) or a Visual studio solution (Windows).
Other options are available.

This should work on most platforms: Linux Windows and Mac.

I have yet to test this on a Mac, so tell me if it works.

There are a number of examples. The simplest is basic_mesh

```
$ examples/basic_mesh/basic_mesh 
building a cube
writing cube.fbx
building a sphere
writing sphere.fbx
building a cylinder
writing cylinder.fbx
building a composite object
writing composite.fbx
```

This writes a number of FBX files for primitive geometry.

The molecules example builds marching cubes or primitive meshes for the Bioblox project.

```
$ examples/molecules/molecules se ../examples/data/2PTC.pdb 
chains EI
60 x 73 x 54
building solvent acessible mesh by inflating the atoms
[0 0][1 5][2 4]...
building solvent excluded mesh by deflating the acessible mesh
[0 0][2 3][1 2]...
writing 2PTC_EI_se_1.0.fbx (13040 vertices)
```


