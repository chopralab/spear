include(ExternalProject)

ExternalProject_Add( CHEMFILES
    GIT_REPOSITORY https://github.com/chemfiles/chemfiles.git
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/chemfiles_build
    CMAKE_CACHE_ARGS    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/chemfiles
                        -DBUILD_SHARED_LIBS:BOOL=ON
                        -DCMAKE_BUILD_TYPE:STRING=Release
)

add_library(chemfiles SHARED IMPORTED)
set_property(TARGET chemfiles PROPERTY IMPORTED_LOCATION
    ${CMAKE_CURRENT_BINARY_DIR}/chemfiles/lib/${CMAKE_SHARED_LIBRARY_PREFIX}chemfiles${CMAKE_SHARED_LIBRARY_SUFFIX}
)
set_property(TARGET chemfiles PROPERTY IMPORTED_IMPLIB
    ${CMAKE_CURRENT_BINARY_DIR}/chemfiles/lib/chemfiles.lib
)

add_dependencies(chemfiles CHEMFILES)

set(CHEMFILES_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/chemfiles/include)
set(CHEMFILES_LIBRARY chemfiles)
