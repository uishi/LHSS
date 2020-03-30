cmake_minimum_required(VERSION 3.10.0)

project(LHSS VERSION 0.0.1 LANGUAGES CXX C)

# TODO
# Make Release mode by default

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(LINK_LIBS gmpxx gmp)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++17 -maes -mbmi2 -O2")

add_subdirectory(demo)
