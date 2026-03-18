cmake_minimum_required(VERSION 3.10)
project(liblombscargle LANGUAGES C CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(lombscargle SHARED
    liblombscargle.cc
    ComplexLS.c
)

target_include_directories(lombscargle PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(NOT WIN32)
    target_link_libraries(lombscargle PRIVATE m)
endif()

if(UNIX AND NOT APPLE)
    # Catch unresolved symbols during build rather than at import-time dlopen.
    target_link_options(lombscargle PRIVATE "-Wl,--no-undefined")
endif()

set_target_properties(lombscargle PROPERTIES
    OUTPUT_NAME "lombscargle"
)

if(WIN32)
    set_target_properties(lombscargle PROPERTIES
        PREFIX "lib"
        WINDOWS_EXPORT_ALL_SYMBOLS ON
    )
endif()
