set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_try test_try.cpp)

enable_testing()
find_package(GTest CONFIG REQUIRED)
target_link_libraries(test_try PRIVATE GTest::gtest GTest::gtest_main GTest::gmock GTest::gmock_main)

add_test(
  NAME google_test
  COMMAND $<TARGET_FILE:test_try>
  )
