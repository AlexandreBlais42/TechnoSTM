cmake_minimum_required(VERSION 3.0)
project(TechnoSTM)

# Version standard C++
set(CMAKE_CXX_STANDARD 17)

add_compile_options(-Wall)

# Définition de variables
set(SOURCE_DIR src/)
set(INCLUDE_DIR src/include/)

# Tous les fichiers sources
file(GLOB SOURCES ${SOURCE_DIR}/*.cpp)

# Répertoire des headers
include_directories(${INCLUDE_DIR})

# Compilation
add_executable(${PROJECT_NAME} ${SOURCES})

