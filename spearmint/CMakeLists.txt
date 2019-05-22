add_library(spearmint SHARED spearmint.cpp read_write.cpp ligand.cpp receptor.cpp)
set_target_properties(spearmint PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(spearmint PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include)

target_compile_definitions(spearmint PRIVATE SPEARMINT_EXPORTS)
set_target_properties(spearmint PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

include(GenerateExportHeader)
generate_export_header(spearmint
    BASE_NAME SPEARMINT
    EXPORT_FILE_NAME ${PROJECT_BINARY_DIR}/include/spear/spearmint_exports.hpp
)

if (MSVC)
	target_link_libraries(spearmint spear chemfiles ws2_32)
else()
	target_link_libraries(spearmint spear chemfiles)
endif()

add_executable(read_and_score read_and_score.c)
target_link_libraries(read_and_score spearmint spear chemfiles)
target_include_directories(read_and_score PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include)

add_executable(read_write_complex read_write_complex.c)
target_link_libraries(read_write_complex spearmint spear chemfiles)
target_include_directories(read_write_complex PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include)

add_executable(print_coords print_coords.c)
target_link_libraries(print_coords spearmint spear chemfiles)
target_include_directories(print_coords PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include)

add_executable(print_bonds print_bonds.c)
target_link_libraries(print_bonds spearmint spear chemfiles)
target_include_directories(print_bonds PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include)

add_executable(print_neighbor print_neighbor.c)
target_link_libraries(print_neighbor spearmint spear chemfiles)
target_include_directories(print_neighbor PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include)
