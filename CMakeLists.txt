cmake_minimum_required(VERSION 3.12)
project(Tadah.MODELS)

# Guard against in-source builds
# Also, make sure user does not accidently invoke cmake .. from the source dir
# of the module in case module is added as a submodule.
get_filename_component(MODULE_TEST ${CMAKE_BINARY_DIR} NAME)
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR}
        OR ${MODULE_TEST} STREQUAL "MODELS"
        )
    message(FATAL_ERROR "
    In-source builds not allowed.
    Please make a new directory (called a build directory) and run CMake from there.
    You may need to remove CMakeCache.txt and CMakeFiles dir. ")
endif()

# check required modules are present
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../CORE")
    message("This module requires CORE module for compilation.")
    message(FATAL_ERROR "CORE module not found.")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

file(GLOB MODELS_SRC CONFIGURE_DEPENDS "*.cpp" "*.h" "*.hpp")
file(GLOB_RECURSE CUTOFFS_SRC CONFIGURE_DEPENDS 
    "cutoffs/*.cpp" "cutoffs/*.h" "cutoffs/*.hpp")
file(GLOB_RECURSE DESCRIPTORS_SRC CONFIGURE_DEPENDS 
    "descriptors/*.cpp" "descriptors/*.h" "descriptors/*.hpp")
file(GLOB_RECURSE FUNCTION_SRC CONFIGURE_DEPENDS 
    "functions/*.cpp" "functions/*.h" "functions/*.hpp")

add_library(tadah.models.obj OBJECT
    ${CUTOFFS_SRC}
    ${FUNCTION_SRC}
    ${DESCRIPTORS_SRC}
    ${MODELS_SRC}
    )

if(TADAH_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(tadah.models.obj PUBLIC OpenMP::OpenMP_CXX)
endif()

add_library(tadah.models STATIC $<TARGET_OBJECTS:tadah.models.obj>)

include(CTest) 
if(TADAH_BUILD_TESTS)
    add_subdirectory(tests)
endif()

set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests/catch2)
