cmake_minimum_required (VERSION 2.6)

project(gilgamesh)


IF(WIN32)
ELSE(WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ENDIF(WIN32)


include_directories(${PROJECT_SOURCE_DIR}/include/ ${PROJECT_SOURCE_DIR}/external/andyzip/include/ ${PROJECT_SOURCE_DIR}/external/glm/)

add_subdirectory(test)
add_subdirectory(examples)

