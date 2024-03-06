Include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.13.7
    )

FetchContent_MakeAvailable(Catch2)

file(GLOB tests_src "*.cpp")

add_executable(test_models ${tests_src})
target_link_libraries(test_models PRIVATE Catch2::Catch2 tadah.models tadah.core)

include(Catch.cmake)
catch_discover_tests(test_models)
