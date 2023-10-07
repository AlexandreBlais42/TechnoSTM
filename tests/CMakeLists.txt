# Compilation
file(GLOB _PROGRAMS *.cpp)
foreach(_PROGRAM ${_PROGRAMS})
  get_filename_component(PROGRAM_NAME ${_PROGRAM} NAME_WE)
  add_executable(${PROGRAM_NAME} ${_PROGRAM})

  target_link_libraries(${PROGRAM_NAME} LINK_PUBLIC STM)
  if (${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL armv7l)
    target_link_libraries(${PROGRAM_NAME} LINK_PUBLIC wiringPi)
  endif()
endforeach()

