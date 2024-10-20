cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_FLAGS "-std=c++11")
project(pro)

include_directories(
    ${PROJECT_SOURCE_DIR}/src/particle
    ${PROJECT_SOURCE_DIR}/src/import
)

set(SRC
    ${PROJECT_SOURCE_DIR}/src/particle.cc
    ${PROJECT_SOURCE_DIR}/src/import.cc
    ${PROJECT_SOURCE_DIR}/src/main.cc
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

ADD_EXECUTABLE(prop ${SRC})