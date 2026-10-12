# Catch2 and its test discovery come from cmake-tools (tadah_setup_tests()).
tadah_setup_tests()

file(GLOB tests_src CONFIGURE_DEPENDS "*.cpp")
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests_data DESTINATION .)

add_executable(test_models ${tests_src})
# Threads: test_kernel_settings_once.cpp builds kernels in std::threads.
find_package(Threads REQUIRED)
target_link_libraries(test_models PRIVATE Catch2 tadah.models tadah.core
                                          Threads::Threads)
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

# tadah.models_ream_goldens_native_link. Pins: the DM_REAM goldens (the
# [ream][rhoform][goldens] cases of test_ream_goldens.cpp) pass in a program
# compiled and linked with -march=native, under LTO when the build has it
# (CMAKE_INTERPROCEDURAL_OPTIMIZATION), as the tadah executable is built
# with TADAH_MARCH_NATIVE: MODELS' code, compiled with MODELS' flags, gives
# such a program the density of every tadah before 1.6.0 bit for bit. GCC
# generates a function's code for the target it was compiled for, at an LTO
# link with -march=native too; a density compiled into the program itself (a
# ReamDensity method moved into a header) fuses its multiply-adds into FMAs
# on an FMA machine, and the goldens fail. Not built when the compiler has
# no -march=native (TADAH_MARCH_NATIVE is ignored then too).
# Before the fix: passes as well (a pin of the build layout).
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-march=native TADAH_MODELS_TESTS_HAVE_MARCH_NATIVE)
if(TADAH_MODELS_TESTS_HAVE_MARCH_NATIVE)
  add_executable(test_models_ream_native_link
      test_ream_goldens.cpp test_main.cpp)
  target_link_libraries(test_models_ream_native_link
      PRIVATE Catch2 tadah.models tadah.core)
  target_compile_options(test_models_ream_native_link PRIVATE -march=native)
  target_link_options(test_models_ream_native_link PRIVATE -march=native)
  add_test(NAME tadah.models_ream_goldens_native_link
      COMMAND test_models_ream_native_link "[ream][rhoform][goldens]"
              --warn NoTests)
else()
  message(STATUS "Tadah.MODELS tests: the compiler has no -march=native; "
                 "tadah.models_ream_goldens_native_link is not built")
endif()

# tadah.models_headers (ctest -L headers): every header of include/
# compiles on its own and twice in a row (cmake-tools'
# tadah_add_header_check). Before the fix: fails on
# descriptors/dm/dm_mEAD.hpp and dm_poly.hpp, which have no include guard
# ('redefinition of ...').
tadah_add_header_check(tadah.models)
