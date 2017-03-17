cmake_minimum_required(VERSION 3.6)
project(SM)

SET(CMAKE_CXX_STANDARD 11)
SET(SOURCE_FILES main.cpp basic.cpp basic.h analytic.cpp analytic.h heston.cpp heston.h plot.cpp plot.h)
SET(LIBRARIES libboost_system.a libboost_iostreams.a libboost_filesystem.a)

INCLUDE_DIRECTORIES(/usr/include/gnuplot-iostream)
LINK_DIRECTORIES(usr/local/lib)

ADD_EXECUTABLE(SM ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(SM ${LIBRARIES})