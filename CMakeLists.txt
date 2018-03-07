cmake_minimum_required(VERSION 3.5)

PROJECT(proj5)
SET(VTK_DIR ~/CIS410SciVis/build-VTK)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj7 MACOSX_BUNDLE proj7)

SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(proj7 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj7 vtkHybrid)
endif()


