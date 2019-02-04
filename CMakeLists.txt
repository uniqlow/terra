cmake_minimum_required(VERSION 3.4 FATAL_ERROR)

project("terra" CXX)
set(CMAKE_CXX_STANDARD 11)

option(TERRA_TEST "Build tests" OFF)

include_directories(${CMAKE_SOURCE_DIR}/include)

if(TERRA_TEST)
	enable_testing()
	add_subdirectory(test)
endif(TERRA_TEST)
