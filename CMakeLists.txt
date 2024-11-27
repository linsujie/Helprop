cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_FLAGS "-std=c++20 -O2")
project(HelProp)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
FIND_PACKAGE(NLOPT)
INCLUDE_DIRECTORIES(${NLOPT_INCLUDE_DIRS})

ADD_SUBDIRECTORY(mongo-c-driver)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/mongo-c-driver/src/libbson/src/bson)

ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/docopt)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/docopt)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${NLOPT_INCLUDE_DIR})

ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/reflect-cpp)
#SET(REFLECTCPP_LIBRARIES ${PROJECT_SOURCE_DIR}/libreflectcpp.a)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/reflect-cpp/include)

IF (PROFILE_DEBUG)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg")
ENDIF()

ADD_LIBRARY(helprop SHARED src/particle.cc src/IO.cc src/fcache.cc src/newton_solver.cc src/root_finding.cc)
TARGET_LINK_LIBRARIES(helprop ${NLOPT_LIBRARIES} ${LIBBSON_LIBRARIES})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

FOREACH(var HelProp iotest)
  ADD_EXECUTABLE(${var} src/${var}.cc)
  TARGET_LINK_LIBRARIES(${var} docopt helprop reflectcpp)
ENDFOREACH()

ADD_EXECUTABLE(get_HCS_distance_test src/get_HCS_distance_test.cc src/particle.cc src/IO.cc src/fcache.cc src/newton_solver.cc src/root_finding.cc)
TARGET_LINK_LIBRARIES(get_HCS_distance_test docopt ${NLOPT_LIBRARIES}
  ${LIBBSON_LIBRARIES} reflectcpp)
