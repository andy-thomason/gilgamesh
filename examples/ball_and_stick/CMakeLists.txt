cmake_minimum_required (VERSION 2.6)

project (ball_and_stick)

add_executable(ball_and_stick main.cpp)

add_definitions(-DCMAKE_SOURCE="${CMAKE_SOURCE_DIR}")
#add_definitions(-DCMAKE_SOURCE="xyz")

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(ball_and_stick Threads::Threads)


