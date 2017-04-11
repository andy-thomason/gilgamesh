cmake_minimum_required (VERSION 2.6)

project (gilgamol)

add_executable(gilgamol main.cpp)

if (WIN32)
else()
  install(TARGETS gilgamol DESTINATION bin)
endif()

add_definitions(-DCMAKE_SOURCE="${CMAKE_SOURCE_DIR}")

find_package(Threads REQUIRED)
find_package (Threads)
target_link_libraries (gilgamol ${CMAKE_THREAD_LIBS_INIT})

if (WIN32)
else()
  set(CPACK_PACKAGE_VERSION "0.9.0")
  set(CPACK_GENERATOR "RPM;DEB")
  set(CPACK_PACKAGE_NAME "gilgamol")
  set(CPACK_PACKAGE_RELEASE 1)
  set(CPACK_PACKAGE_CONTACT "Andy Thomason")
  set(CPACK_PACKAGE_VENDOR "Atomic Increment")
  set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_PACKAGE_RELEASE}.${CMAKE_SYSTEM_PROCESSOR}")
  SET(CPACK_PACKAGE_EXECUTABLES "gilgamol" "gilgamol")
  include(CPack)
endif()
