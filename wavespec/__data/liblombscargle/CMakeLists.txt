cmake_minimum_required(VERSION 3.10)
project(liblombscargle LANGUAGES C CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(lombscargle SHARED
    liblombscargle.cc
    ComplexLS.c
)

target_include_directories(lombscargle PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(lombscargle PRIVATE m)

set_target_properties(lombscargle PROPERTIES
    OUTPUT_NAME "lombscargle"
)
