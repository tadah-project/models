# Catch2 and its test discovery come from cmake-tools (tadah_setup_tests()).
tadah_setup_tests()

file(GLOB tests_src CONFIGURE_DEPENDS "*.cpp")
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests_data DESTINATION .)

add_executable(test_models ${tests_src})
target_link_libraries(test_models PRIVATE Catch2 tadah.models tadah.core)
target_include_directories(test_models PRIVATE ../include)

include(Catch)
# WORKING_DIRECTORY pins each registered Catch test case's cwd to the
# binary's directory so tests opening fixtures via relative paths
# (e.g. "tests_data/basis", "tests_data/eam.alloy" in test_kernels.cpp /
# test_dm.cpp / test_d2.cpp) resolve under ctest and IDE test runners
# regardless of the caller's cwd. Direct ./test_models invocations from
# elsewhere still need to cd into the binary dir.
catch_discover_tests(test_models
    WORKING_DIRECTORY $<TARGET_FILE_DIR:test_models>
)
