function(spear_define_test _name_)
    add_test(NAME ${_name_}
        COMMAND ${RUNNER_COMMAND} ${_name_}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    if(WIN32)
        # We need to set the path to allow access to spear.dll
        # (and any other DLL)
        STRING(REPLACE ";" "\\;" PATH_STRING "$ENV{PATH}")
        set_tests_properties(${_name_}
            PROPERTIES ENVIRONMENT "PATH=${PATH_STRING}\;$<TARGET_FILE_DIR:spear>\;${CMAKE_BINARY_DIR}\\chemfiles\\bin"
        )
    endif()
endfunction()

function(spear_unit_test _name_ _file_)
    add_executable(${_name_} ${_file_} ${SPEAR_OBJECTS})
    set_target_properties(${_name_} PROPERTIES LINKER_LANGUAGE CXX)
    target_include_directories(${_name_} SYSTEM PRIVATE
        ${CMAKE_SOURCE_DIR}/tests/external
    )

    spear_define_test(${_name_})
endfunction()

function(spear_cpp_test _file_)
    get_filename_component(_name_ ${_file_} NAME_WE)
    spear_unit_test(${_name_} ${_file_})
	if (MSVC)
		target_link_libraries(${_name_} spear chemfiles ws2_32)
	else()
		target_link_libraries(${_name_} spear chemfiles)
	endif()
endfunction()

if (NOT MSVC)
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/share share
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	)

	execute_process(
		COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data data
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	)
else()
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/share share
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	)

	execute_process(
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data data
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	)
endif()

file(GLOB all_test_files
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

foreach(test_file IN LISTS all_test_files)
    spear_cpp_test(${test_file})
endforeach(test_file)
