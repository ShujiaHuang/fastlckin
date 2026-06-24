# ── tests/CMakeLists.txt ────────────────────────────────────────────
# Unit tests and test data generator for fastlckin
# This file is included from the top-level CMakeLists.txt via add_subdirectory(tests)

# Variables inherited from parent scope:
#   HTSLIB_A, HTSLIB_DIR, COMMON_SYS_LIB, Buildhts (target)
#   CMAKE_CURRENT_SOURCE_DIR here = <project>/tests

set(PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)

# ── Shared library sources (everything except main.cpp) ─────────────
file(GLOB_RECURSE LIB_SRCS "${PROJECT_ROOT}/src/*.cpp")
list(FILTER LIB_SRCS EXCLUDE REGEX ".*main\\.cpp$")

# ── Shared static library: compile once, link from all tests ────────
# This avoids each test target re-compiling the full source tree,
# dramatically reducing CI memory usage during parallel builds.
add_library(fastlckin_lib STATIC ${LIB_SRCS})
add_dependencies(fastlckin_lib Buildhts)

target_include_directories(fastlckin_lib PUBLIC
    ${PROJECT_ROOT}/src
    ${HTSLIB_DIR}
    ${HTSLIB_DIR}/htslib
)

if(APPLE)
    target_link_options(fastlckin_lib PUBLIC -Wl,-no_compact_unwind)
endif()

target_link_libraries(fastlckin_lib ${HTSLIB_A} ${COMMON_SYS_LIB})

# ── Test data generator (standalone, no htslib dependency) ──────────
add_executable(generate_test_data ${CMAKE_CURRENT_SOURCE_DIR}/generate_test_data.cpp)
set_target_properties(generate_test_data PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests
)
add_test(NAME generate_test_data
         COMMAND ${CMAKE_BINARY_DIR}/tests/generate_test_data
         WORKING_DIRECTORY ${PROJECT_ROOT})

# ── Known-relationship data generator (standalone, no htslib) ──────
add_executable(generate_known_relationships
    ${CMAKE_CURRENT_SOURCE_DIR}/generate_known_relationships.cpp)
set_target_properties(generate_known_relationships PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests
)
add_test(NAME generate_known_relationships
         COMMAND ${CMAKE_BINARY_DIR}/tests/generate_known_relationships
         WORKING_DIRECTORY ${PROJECT_ROOT})

# ── Unit test executables ───────────────────────────────────────────
# Each test only compiles its own .cpp and links against fastlckin_lib
file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/test_*.cpp")

foreach(TEST_SRC ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_SRC} NAME_WE)

    add_executable(${TEST_NAME} ${TEST_SRC})
    target_link_libraries(${TEST_NAME} fastlckin_lib)

    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME}
             WORKING_DIRECTORY ${PROJECT_ROOT})
endforeach()
