Include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.13.7
    )

FetchContent_MakeAvailable(Catch2)

file(GLOB tests_src "*.cpp")

find_library(tadah.core NAMES tadah.core
    PATHS ${CMAKE_CURRENT_SOURCE_DIR}/../../CORE/*
    )
if(tadah.core)
    message(STATUS "Found ${tadah.core}") 
else()
    message(FATAL_ERROR "
    Could not find ltadah.core in the CORE module. 
    Hint: Recompile CORE.
    ")
endif()

add_executable(test_models ${tests_src})
target_link_libraries(test_models PRIVATE Catch2::Catch2 tadah.models ${tadah.core})

include(Catch.cmake)
catch_discover_tests(test_models)
