add_library(spearmint spearmint.cpp read_write.cpp)
set_target_properties(spearmint PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(spearmint PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include)

if (MSVC)
	target_link_libraries(spearmint spear chemfiles ws2_32)
else()
	target_link_libraries(spearmint spear chemfiles)
endif()

message(STATUS ${PROJECT_BINARY_DIR}/include)

add_executable(read_and_score read_and_score.c)
target_link_libraries(read_and_score spearmint spear chemfiles)
target_include_directories(read_and_score PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include)

add_executable(read_write_complex read_write_complex.c)
target_link_libraries(read_write_complex spearmint spear chemfiles)
target_include_directories(read_write_complex PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include)
