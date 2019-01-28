cmake_minimum_required(VERSION 3.8)
project(p32)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

include_directories(${CMAKE_SOURCE_DIR}/src)

set(SOURCE_FILES
        ${CMAKE_SOURCE_DIR}/src/input_image.h
        ${CMAKE_SOURCE_DIR}/src/input_image.cc
        ${CMAKE_SOURCE_DIR}/src/complex.h
        ${CMAKE_SOURCE_DIR}/src/complex.cc
        )

add_executable(${PROJECT_NAME} ${SOURCE_FILES} main_MPI.cc)
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})
