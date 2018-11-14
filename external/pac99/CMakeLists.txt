project(PAC99 Fortran)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")

add_executable(pac99 ${PROJECT_SOURCE_DIR}/pac99.f)

install(TARGETS pac99 DESTINATION bin)
