execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_SOURCE_DIR}/csd_distributions.dat.tgz
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
