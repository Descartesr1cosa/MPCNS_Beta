cmake_minimum_required(VERSION 3.5)
project(MPCNS_HallMHD)
set(CMAKE_CXX_STANDARD 20)


# 从环境变量里读 PETSC_DIR / PETSC_ARCH
if(NOT DEFINED ENV{PETSC_DIR})
    message(FATAL_ERROR "环境变量 PETSC_DIR 没设，例如 export PETSC_DIR=/home/Descartes/tmp_petsc/petsc")
endif()

if(NOT DEFINED ENV{PETSC_ARCH})
    message(FATAL_ERROR "环境变量 PETSC_ARCH 没设，例如 export PETSC_ARCH=arch-mswin-c-debug 或 arch-mswin-c-opt")
endif()

set(PETSC_DIR  $ENV{PETSC_DIR})
set(PETSC_ARCH $ENV{PETSC_ARCH})

message(STATUS "Using PETSc from: ${PETSC_DIR}")
message(STATUS "Using PETSc arch: ${PETSC_ARCH}")

# PETSc 头文件
include_directories(
    "${PETSC_DIR}/include"
    "${PETSC_DIR}/${PETSC_ARCH}/include"
)

# PETSc 库目录
link_directories("${PETSC_DIR}/${PETSC_ARCH}/lib")


# Find MPI
find_package(MPI REQUIRED)

# Set Compile Options
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# include files directory
include_directories(${MPI_INCLUDE_PATH})
include_directories(../include)
include_directories(../Z3_HallMHD)
include_directories(.)

# source files directiory
aux_source_directory(../src/0_basic  DIRSRCS0)
aux_source_directory(../src/1_grid  DIRSRCS1)
aux_source_directory(../src/2_topology  DIRSRCS2)
aux_source_directory(../src/3_field  DIRSRCS3)

aux_source_directory(./  DIRSRCS4)

# build the executable program
add_executable(MPCNS_HallMHD main.cpp ${DIRSRCS0} ${DIRSRCS1} ${DIRSRCS2} ${DIRSRCS3} ${DIRSRCS4})

# set compile options
#target_compile_options(MPCNS_HallMHD PRIVATE -std=c++20 -fconcepts)

# add the mpi linkings
target_link_libraries(MPCNS_HallMHD PRIVATE
    petsc     openblas
${MPI_LIBRARIES} stdc++fs)
