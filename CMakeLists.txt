cmake_minimum_required(VERSION 2.8.12)

project(spear C CXX)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake")
if (POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW) # OSX RPATH handling
endif()

file(READ ${PROJECT_SOURCE_DIR}/VERSION SPEAR_VERSION)
string(STRIP ${SPEAR_VERSION} SPEAR_VERSION)
string(REGEX REPLACE "^([0-9]+)\\..*" "\\1" SPEAR_VERSION_MAJOR "${SPEAR_VERSION}")
string(REGEX REPLACE "^[0-9]+\\.([0-9]+).*" "\\1" SPEAR_VERSION_MINOR "${SPEAR_VERSION}")
string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" SPEAR_VERSION_PATCH "${SPEAR_VERSION}")
set(SPEAR_VERSION_SHORT "${SPEAR_VERSION_MAJOR}.${SPEAR_VERSION_MINOR}.${SPEAR_VERSION_PATCH}")

option(BUILD_SHARED_LIBS "Build shared libraries instead of static ones" OFF)

option(SPEAR_BUILD_TESTS "Build unit tests." OFF)
option(SPEAR_BUILD_DOCUMENTATION "Build the documentation." OFF)
option(SPEAR_CODE_COVERAGE "Enable code coverage" OFF)
option(SPEAR_USE_WARNINGS "Compile the code with warnings (default in debug mode)" OFF)
option(SPEAR_BUILD_CHEMFILES "Build chemfiles in addition to spear." ON)

mark_as_advanced(SPEAR_BUILD_DOCTESTS)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'release' as none was specified.")
    set(CMAKE_BUILD_TYPE "release"
        CACHE STRING
        "Choose the type of build, options are: none(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) debug release relwithdebinfo minsizerel."
    FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS release debug relwithdebuginfo minsizerel none)
endif()

include(CompilerFlags)

if(SPEAR_USE_WARNINGS)
    # Make sure to have asserts compiled
    add_definitions("-UNDEBUG")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SPEAR_CXX_WARNINGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SPEAR_C_WARNINGS}")
endif()
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${SPEAR_CXX_WARNINGS}")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${SPEAR_C_WARNINGS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SPEAR_SANITIZERS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SPEAR_SANITIZERS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${SPEAR_SANITIZERS}")

if(WIN32)
    add_definitions("-DSPEAR_WINDOWS")
endif()

if(MSVC)
    add_definitions("-D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_WARNINGS")
endif()

if(SPEAR_CODE_COVERAGE)
    message(STATUS "Code coverage enabled")
    if(NOT CMAKE_COMPILER_IS_GNUCXX)
        message(FATAL_ERROR "Code coverage can only be used with GCC")
    endif()
    # Code coverage should use gcc
    set(CMAKE_BUILD_TYPE "debug" CACHE STRING "Type of build: release|debug|relwithdebuginfo" FORCE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
endif()

# Attempt to find chemfiles
if(NOT SPEAR_BUILD_CHEMFILES)
    find_package(chemfiles REQUIRED)
else()
    message(STATUS "Chemfiles is being built!")
    include(InstallChemfiles)
endif()

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
find_package(Boost REQUIRED COMPONENTS graph filesystem program_options)
find_package(Eigen3 REQUIRED)

file(GLOB_RECURSE SPEAR_objects_sources src/**.cpp)
add_library(spear ${SPEAR_objects_sources})
target_include_directories(spear PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(spear PUBLIC
    Boost::graph
    chemfiles
)

target_compile_definitions(spear PRIVATE SPEAR_EXPORTS)
if(${BUILD_SHARED_LIBS})
    set_target_properties(spear PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
endif()

set_property(TARGET spear PROPERTY VERSION ${SPEAR_VERSION_SHORT})
set_property(TARGET spear PROPERTY SOVERSION ${SPEAR_VERSION_SHORT})

configure_file (
  "${PROJECT_SOURCE_DIR}/include/spear/config.in.hpp"
  "${PROJECT_BINARY_DIR}/include/spear/config.hpp"
)

include(GenerateExportHeader)
generate_export_header(spear
    BASE_NAME SPEAR
    EXPORT_FILE_NAME ${PROJECT_BINARY_DIR}/include/spear/exports.hpp
)

#----------------------------------------------------------------------------------------#
# Installation configuration
#----------------------------------------------------------------------------------------#
set(INCLUDE_INSTALL_DIR "include")
install(TARGETS spear EXPORT spear-targets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib${LIB_SUFFIX}
    ARCHIVE DESTINATION lib${LIB_SUFFIX}
    INCLUDES DESTINATION include
)
install(DIRECTORY include/ DESTINATION include)
install(DIRECTORY ${PROJECT_BINARY_DIR}/include/ DESTINATION include)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/spear-config-version.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/spear-config-version.cmake"
    @ONLY
)

# Install package files
set(PACKAGE_LOCATION lib${LIB_SUFFIX}/cmake/spear)
install(EXPORT spear-targets
  FILE spear-targets.cmake
  DESTINATION ${PACKAGE_LOCATION}
)
install(FILES
    cmake/spear-config.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/spear-config-version.cmake"
    DESTINATION ${PACKAGE_LOCATION}
    COMPONENT Devel
)

# uninstall target
configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/uninstall.in.cmake"
    "${PROJECT_BINARY_DIR}/uninstall.cmake"
    @ONLY
)
add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${PROJECT_BINARY_DIR}/uninstall.cmake
)

add_subdirectory(share)

if(SPEAR_BUILD_DOCUMENTATION)
    add_subdirectory(doc)
endif()

enable_testing()
if(SPEAR_BUILD_TESTS)
    add_subdirectory(tests)
endif()
