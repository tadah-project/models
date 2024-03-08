file(GLOB tests_src "*.cpp")

add_executable(test_models ${tests_src})
target_link_libraries(test_models PRIVATE Catch2 tadah.models tadah.core)

include(Catch.cmake)
catch_discover_tests(test_models)
