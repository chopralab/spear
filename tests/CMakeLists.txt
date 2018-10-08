function(spear_define_test _name_)
    add_test(NAME ${_name_}
        COMMAND ${RUNNER_COMMAND} ./${_name_}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    if(WIN32)
        # We need to set the path to allow access to spear.dll
        # (and any other DLL)
        STRING(REPLACE ";" "\\;" PATH_STRING "$ENV{PATH}")
        set_tests_properties(${_name_}
            PROPERTIES ENVIRONMENT "PATH=${PATH_STRING}\;$<TARGET_FILE_DIR:spear>"
        )
    endif()
endfunction()

function(spear_unit_test _name_ _file_)
    add_executable(${_name_} ${_file_} ${SPEAR_OBJECTS})
    # We need to pretend we are inside the DLL to access all functions
    target_compile_definitions(${_name_} PRIVATE spear_EXPORTS)
    set_target_properties(${_name_} PROPERTIES LINKER_LANGUAGE CXX)
    target_include_directories(${_name_} SYSTEM PRIVATE
        ${BOOST_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_BINARY_DIR}/include
        ${CMAKE_BINARY_DIR}/chemfiles/include
        ${CMAKE_SOURCE_DIR}/tests/external
    )

    spear_define_test(${_name_})
endfunction()

function(spear_cpp_test _file_)
    get_filename_component(_name_ ${_file_} NAME_WE)
    spear_unit_test(${_name_} ${_file_})
    target_link_libraries(${_name_} spear chemfiles)
endfunction()

file(GLOB all_test_files
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

foreach(test_file IN LISTS all_test_files)
    spear_cpp_test(${test_file})
endforeach(test_file)