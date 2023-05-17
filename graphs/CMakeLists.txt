cmake_minimum_required (VERSION 3.0)

set(CMAKE_CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED ON)

project (GRAPH)

set(OpenGL_GL_PREFERENCE "GLVND")
find_package( OpenGL REQUIRED )
find_package( glm CONFIG REQUIRED )
find_package( SDL2 CONFIG REQUIRED )

include_directories(SYSTEM anttweakbar/include)

add_subdirectory(anttweakbar EXCLUDE_FROM_ALL)

set(target graph.elf)

set(graph_SOURCES
  shader_utils.cpp
  libUtils.cpp
  libFFTcstm.cpp
  libCZTcstm.cpp
  fffc.cpp
  graph.cpp
)

add_executable(${target} ${graph_SOURCES})

target_link_libraries(${target}
  PRIVATE
  GL
  SDL2
  AntTweakBar
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/graph.f.glsl DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/graph.v.glsl DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
