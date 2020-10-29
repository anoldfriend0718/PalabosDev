set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set(TEST_ASSEMBLY "tests")
file(GLOB_RECURSE TEST_CPPS ${CUSTOMIZED_PARABOS_SOURCE_REPO}/tests/main.cpp ${CUSTOMIZED_PARABOS_SOURCE_REPO}/tests/test_*.cpp)
add_executable(${TEST_ASSEMBLY} ${TEST_CPPS})


find_package(GTest CONFIG REQUIRED)
# Link with the following libraries
target_link_libraries(${TEST_ASSEMBLY} PRIVATE GTest::gtest GTest::gtest_main GTest::gmock GTest::gmock_main)
target_link_libraries(${TEST_ASSEMBLY} PRIVATE ${PALABOS_LIB_NAME})

if(ENABLE_MPI)
    target_link_libraries(${TEST_ASSEMBLY} PUBLIC ${MPI_CXX_LIBRARIES})
endif()

add_test(
  NAME google_test
  COMMAND $<TARGET_FILE:${TEST_ASSEMBLY}>
  )
