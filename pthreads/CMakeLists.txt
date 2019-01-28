cmake_minimum_required(VERSION 3.9)
project(p31)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -O3 -pthread")

include_directories(
    ${CMAKE_SOURCE_DIR}/src
)

set(SOURCE_FILES
    ${CMAKE_SOURCE_DIR}/src/complex.h
    ${CMAKE_SOURCE_DIR}/src/complex.cc
)


add_executable(${PROJECT_NAME} ${SOURCE_FILES} main.cc)
