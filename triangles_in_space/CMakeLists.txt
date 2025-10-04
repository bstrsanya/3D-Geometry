cmake_minimum_required(VERSION 3.10)

project(triangle_calculations)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")
set(CMAKE_CXX_FLAGS "-O3")
add_subdirectory(tests)
include_directories(${CMAKE_SOURCE_DIR}/include)
add_executable(${PROJECT_NAME} main.cpp)
# target_link_libraries(${PROJECT_NAME} triangle)