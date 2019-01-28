cmake_minimum_required(VERSION 3.9)
project(p33)

find_package(CUDA QUIET REQUIRED)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-std=c++11;)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS})

include_directories(
${CMAKE_SOURCE_DIR}/src
)
set(SOURCE_FILES
${CMAKE_SOURCE_DIR}/src/input_image.h
${CMAKE_SOURCE_DIR}/src/input_image.cc
${CMAKE_SOURCE_DIR}/src/complex.h
${CMAKE_SOURCE_DIR}/src/complex.cc
)
cuda_add_executable(${PROJECT_NAME} ${SOURCE_FILES} main_CUDA.cu)
