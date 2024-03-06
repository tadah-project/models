Include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.13.7
    )

FetchContent_MakeAvailable(Catch2)

file(GLOB tests_src "*.cpp")


include(ExternalProject)
ExternalProject_Add(Tadah.CORE
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../CORE
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../CORE/build
    INSTALL_COMMAND cmake -E echo "Skipping install step."
    )
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../CORE/build)

add_executable(test_models ${tests_src})
target_link_libraries(test_models PRIVATE Catch2::Catch2 tadah.models tadah.core)

include(Catch.cmake)
catch_discover_tests(test_models)
