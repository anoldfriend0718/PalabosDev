
set(EXECUTABLE_NAME "cylinder2d")
set(MODIFIED_EXECUTABLE_NAME "guocylinder2d")
set(TEST_EXECUTABLE_NAME "lyz_guo_cylinder2d")
set(A_EXECUTABLE_NAME "testcylinder2d")
set(B_EXECUTABLE_NAME "poiseuille")

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ###############################################################################
# # Execution
# ###############################################################################

add_executable(${EXECUTABLE_NAME} "./${EXECUTABLE_NAME}.cpp")
add_executable(${MODIFIED_EXECUTABLE_NAME} "./${MODIFIED_EXECUTABLE_NAME}.cpp")
add_executable(${TEST_EXECUTABLE_NAME} "./${TEST_EXECUTABLE_NAME}.cpp")
add_executable(${A_EXECUTABLE_NAME} "./${A_EXECUTABLE_NAME}.cpp")
add_executable(${B_EXECUTABLE_NAME} "./${B_EXECUTABLE_NAME}.cpp")
# Link with the following libraries
target_link_libraries(${EXECUTABLE_NAME} ${PALABOS_LIB_NAME})
target_link_libraries(${MODIFIED_EXECUTABLE_NAME} ${PALABOS_LIB_NAME})
target_link_libraries(${TEST_EXECUTABLE_NAME} ${PALABOS_LIB_NAME})
target_link_libraries(${A_EXECUTABLE_NAME} ${PALABOS_LIB_NAME})
target_link_libraries(${B_EXECUTABLE_NAME} ${PALABOS_LIB_NAME})
if(ENABLE_MPI)
    target_link_libraries(${EXECUTABLE_NAME} ${MPI_CXX_LIBRARIES})
    target_link_libraries(${MODIFIED_EXECUTABLE_NAME} ${MPI_CXX_LIBRARIES})
    target_link_libraries(${TEST_EXECUTABLE_NAME} ${MPI_CXX_LIBRARIES})
    target_link_libraries(${A_EXECUTABLE_NAME} ${MPI_CXX_LIBRARIES})
    target_link_libraries(${B_EXECUTABLE_NAME} ${MPI_CXX_LIBRARIES})
endif()
