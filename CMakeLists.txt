cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(reference-crystal-matching VERSION 1.0.1 LANGUAGES CXX)

find_package(coretoolkit REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(structure-identification REQUIRED)

file(GLOB PTM_MATH_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/dependencies/ptm_math/*.cpp)

volt_add_plugin(reference-crystal-matching
    SOURCES
        src/analysis/reference_lattice.cpp
        src/analysis/grain_frame.cpp
        src/analysis/neighbor_graph.cpp
        src/analysis/site_assignment.cpp
        src/analysis/full_crystal_context.cpp
        src/analysis/reference_crystal_matching_service.cpp
    MAIN_SOURCE src/plugin.cpp
    VENDOR_SOURCES ${PTM_MATH_SOURCES}
    INCLUDE_DIRS dependencies/ptm_math src/analysis
    DEPENDENCIES
        structure-identification::structure-identification
        yaml-cpp::yaml-cpp
    INSTALL_DATA structures
)

# ---- Tests ----------------------------------------------------------------
# Generic suite proving RCM is not forsterite-specific + locking the
# reciprocity invariant. BUILD_TESTING defaults ON so `ctest` works after a
# normal local build; the conan package build() passes -DBUILD_TESTING=OFF.
option(BUILD_TESTING "Build the reference-crystal-matching test suite" ON)
if(BUILD_TESTING)
    enable_testing()
    add_executable(test_rcm tests/test_rcm.cpp)
    # Linking the static lib propagates its PUBLIC includes (include/,
    # dependencies/ptm_math, src/analysis) and all PUBLIC deps. The lib's PCH is
    # PRIVATE, so the test includes <volt/core/volt.h> transitively via headers.
    target_link_libraries(test_rcm PRIVATE reference-crystal-matching_lib TBB::tbb)
    target_compile_features(test_rcm PRIVATE cxx_std_23)
    target_compile_definitions(test_rcm PRIVATE
        RCM_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures")
    add_test(NAME rcm_suite
             COMMAND test_rcm "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures")
endif()
