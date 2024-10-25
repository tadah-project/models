file(GLOB tests_src "*.cpp")
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests_data DESTINATION .)

add_executable(test_models ${tests_src})
target_link_libraries(test_models PRIVATE Catch2 tadah.models tadah.core)
target_include_directories(test_models PRIVATE ../include)

include(Catch.cmake)
catch_discover_tests(test_models)
