cmake_minimum_required(VERSION 2.6)
project(TUR3D)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
set(CMAKE_BUILD_TYPE "Debug")

find_package(HDF5 REQUIRED)
find_package(GSL REQUIRED)

include_directories(${HDF5_INCLUDE_DIR} ${GSL_INCLUDE_DIR} ${CMAKE_SOURCE_DIR})


#add_executable(tur3d.exe tur3d.c types.c io.c)
add_executable(fou3d.exe fou3d.c types.c io.c)

#target_link_libraries(tur3d.exe fftw3 ${GSL_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(fou3d.exe fftw3 ${GSL_LIBRARIES} ${HDF5_LIBRARIES})
