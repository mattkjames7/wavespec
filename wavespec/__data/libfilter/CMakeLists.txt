cmake_minimum_required(VERSION 3.10)
project(libfilter LANGUAGES C CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(filter SHARED
    libfilter.cc
)

target_include_directories(filter PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(filter PRIVATE m)

set_target_properties(filter PROPERTIES
    OUTPUT_NAME "filter"
)

if(WIN32)
    set_target_properties(filter PROPERTIES PREFIX "lib")
endif()
