cmake_minimum_required (VERSION 2.6)

project (molecules)

add_executable(molecules main.cpp)

add_definitions(-DCMAKE_SOURCE="${CMAKE_SOURCE}")

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(molecules Threads::Threads)


