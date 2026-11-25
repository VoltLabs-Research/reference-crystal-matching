cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(reference-crystal-matching VERSION 1.0.1 LANGUAGES CXX)

if(NOT TARGET coretoolkit::coretoolkit)
    find_package(coretoolkit REQUIRED)
endif()
if(NOT TARGET yaml-cpp::yaml-cpp)
    find_package(yaml-cpp REQUIRED)
endif()
if(NOT TARGET structure-identification::structure-identification)
    find_package(structure-identification REQUIRED)
endif()

file(GLOB PTM_MATH_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/ptm_math/*.cpp)

volt_add_plugin(reference-crystal-matching
    SOURCES
        src/registration/reference_lattice.cpp
        src/registration/grain_frame.cpp
        src/registration/grain_segmentation.cpp
        src/registration/neighbor_graph.cpp
        src/registration/site_assignment.cpp
        src/registration/crystal_context.cpp
        src/registration/reference_crystal_matching_service.cpp
        src/registration/reference_crystal_matching_log.cpp
    MAIN_SOURCE src/plugin.cpp
    VENDOR_SOURCES ${PTM_MATH_SOURCES}
    INCLUDE_DIRS dependencies/ptm_math
    DEPENDENCIES
        structure-identification::structure-identification
        yaml-cpp::yaml-cpp
    INSTALL_DATA structures
)
