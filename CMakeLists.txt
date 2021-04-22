cmake_minimum_required(VERSION 3.15)

project(NanoLogLite C CXX)

SET(CMAKE_CXX_COMPILER "g++")
SET(CMAKE_C_COMPILER "gcc")

add_compile_options(-Wall -Wno-sign-compare -Wno-unused-function -Werror -O3 -g -std=c++17)

add_library(nanologlite SHARED NanoLog.cc RuntimeLogger.cc)
install(TARGETS nanologlite)

add_subdirectory(test)
