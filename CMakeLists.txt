cmake_minimum_required (VERSION 3.1.0)

project(gilgamesh)

set(VERSION "0.9.0")
set(CMAKE_CXX_STANDARD 11)

include_directories(${PROJECT_SOURCE_DIR}/include/ ${PROJECT_SOURCE_DIR}/external/andyzip/include/ ${PROJECT_SOURCE_DIR}/external/glm/)

add_subdirectory(test)
add_subdirectory(examples)

