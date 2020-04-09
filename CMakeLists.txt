cmake_minimum_required(VERSION 3.0.0)

project(LHSS VERSION 0.0.1 LANGUAGES CXX C)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
 set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

set(CMAKE_CXX_COMPILERS /usr/bin/g++)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++17 -mbmi2 -O2")

set(LINK_LIBS gmpxx gmp)

add_subdirectory(demo)
