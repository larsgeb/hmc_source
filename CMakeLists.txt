cmake_minimum_required(VERSION 3.7)
project(hmc_source)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        aux.cpp
        aux.hpp
        make_data.cpp
        mc.cpp
        mc.hpp
        sampling.cpp)

add_executable(hmc_source ${SOURCE_FILES})