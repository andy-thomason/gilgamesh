cmake_minimum_required (VERSION 2.6)

project (fbx_dump)

add_definitions(-DCMAKE_SOURCE="${CMAKE_SOURCE_DIR}")

add_executable(fbx_dump main.cpp)
