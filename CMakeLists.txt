cmake_minimum_required(VERSION 3.10)
project(Ray)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


add_executable(Ray src/main.cpp)
target_include_directories(Ray PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)