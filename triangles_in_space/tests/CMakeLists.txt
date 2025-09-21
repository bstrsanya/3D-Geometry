cmake_minimum_required(VERSION 3.10)

project(tests)
add_executable(${PROJECT_NAME} test.cpp)
target_link_libraries (${PROJECT_NAME} PUBLIC
    triangle
    gtest
    gtest_main
)