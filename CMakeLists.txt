cmake_minimum_required(VERSION 3.16)
project(blackhole CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

add_executable(blackhole
  src/main.cpp
  src/render.cpp
  src/spectrum.cpp
  src/image.cpp
  src/validate.cpp
)

target_include_directories(blackhole PRIVATE src)

if(MSVC)
  target_compile_options(blackhole PRIVATE /O2 /fp:precise)
else()
  # -O3 is fine, but -ffast-math is deliberately NOT used: the conservation
  # checks depend on IEEE semantics, and reassociating the geodesic right-hand
  # side would quietly cost accuracy.
  target_compile_options(blackhole PRIVATE -O3 -Wall -Wextra)
endif()

find_package(Threads REQUIRED)
target_link_libraries(blackhole PRIVATE Threads::Threads)
