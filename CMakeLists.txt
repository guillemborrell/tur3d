cmake_minimum_required(VERSION 2.6)
project(TUR3D)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

find_package(HDF5 REQUIRED)
find_package(GSL REQUIRED)
#find_package(FFTW REQUIRED)

include_directories(${HDF5_INCLUDE_DIR} ${GSL_INCLUDE_DIR})

add_executable(tur3d.exe tur3d.c)
target_link_libraries(tur3d.exe fftw3 ${GSL_LIBRARIES} ${HDF5_LIBRARIES})
