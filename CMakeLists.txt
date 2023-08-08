
#cmake_minimum_required (VERSION 3.27)

#For old environments
cmake_minimum_required (VERSION 3.11)


project(pmat LANGUAGES CXX VERSION 0.1 DESCRIPTION "SIMPLE C++ STATIC LIBRARY FOR MATRIX OPERATIONS")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#For old environments
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pthread -Wno-attributes" )

enable_testing()

add_subdirectory (src)
add_subdirectory (test) 




