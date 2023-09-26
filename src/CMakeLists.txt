file(GLOB_RECURSE SOURCES *.cpp)
add_library(STM ${SOURCES})

target_include_directories(STM PUBLIC include/)
