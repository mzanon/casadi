cmake_minimum_required(VERSION 2.8.6)
include_directories(../)

# Linear solver using hsl library
casadi_plugin(LinearSolver hsl
  hsl.hpp
  hsl.cpp
  hsl_meta.cpp)

casadi_plugin_link_libraries(LinearSolver hsl ${HSL_LIBRARIES})

