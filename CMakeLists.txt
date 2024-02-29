cmake_minimum_required(VERSION 3.12)

project(Tadah.MODELS)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB_RECURSE MODELS_SRC CONFIGURE_DEPENDS "*.cpp")

add_library(tadah.models STATIC
    ${MODELS_SRC}
    )

