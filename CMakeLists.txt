cmake_minimum_required(VERSION 3.12)
project(Tadah.MODELS LANGUAGES CXX)

# 0. Guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR
    "In-source builds are not allowed.\n"
    "Create a separate build directory and run CMake from there.")
endif()

# 1. Global settings
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# default build type for single-config generators (Ninja/Make)
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

if(BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)

  # Always full RPATH
  # use, i.e. don't skip the full RPATH for the build tree
  set(CMAKE_SKIP_BUILD_RPATH FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

# 2. Bring in the shared compile/link options
include(FetchContent)
if(NOT TARGET tadah::compile_options)
  FetchContent_Declare (
    tadah_cmake
    GIT_REPOSITORY  https://git.ecdf.ed.ac.uk/tadah/cmake-tools.git
    GIT_TAG         main
  )
  FetchContent_MakeAvailable(tadah_cmake)   # defines tadah::compile_options
endif()

# (option) Developer builds -DTADAH_DEVELOP=ON
tadah_use_local_module(Tadah.CORE CORE)  # FetchContent name + folder name

# 3. Git branch information, for informational purposes 
tadah_get_git_branch(GIT_BRANCH)
message(STATUS "${PROJECT_NAME} module is using GIT_BRANCH: ${GIT_BRANCH}")

#4. Fetch Tadah! modules
FetchContent_Declare(
    Tadah.CORE
    GIT_REPOSITORY https://git.ecdf.ed.ac.uk/tadah/core.git
    GIT_TAG origin/${GIT_BRANCH}
)
FetchContent_MakeAvailable(Tadah.CORE)

file(GLOB MODELS_SRC CONFIGURE_DEPENDS "src/*.cpp")
add_library(tadah.models ${MODELS_SRC})
target_include_directories(tadah.models
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(tadah.models
  PUBLIC
    tadah::compile_options    # get warnings/asan/-pg flags
    tadah.core)

if(TADAH_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(tadah.models PUBLIC OpenMP::OpenMP_CXX)
endif()

set(BUILD_TESTING OFF CACHE BOOL "Build tests")
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()

# Installation rules
install(TARGETS tadah.models
        EXPORT  tadahModelsTargets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(DIRECTORY include/ DESTINATION include)
install(EXPORT tadahModelsTargets
        NAMESPACE tadah::
        DESTINATION lib/cmake/tadah.models)
