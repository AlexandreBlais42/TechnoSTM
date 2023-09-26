cmake_minimum_required(VERSION 3.0)
project(TechnoSTM)

# Version standard C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_compile_options(-Wall -O2)

add_subdirectory(programs)
add_subdirectory(src)
