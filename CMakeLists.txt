cmake_minimum_required(VERSION 3.5)
project(MsaCalc)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

find_package(OpenMP REQUIRED)

add_executable(MsaCalc Test/test.cpp msa.h Gnuplot.h)
#target_compile_options(MsaCalc PRIVATE -g -O3 -Wall -fopenmp)
target_link_libraries(MsaCalc PRIVATE ${OpenMP_CXX_FLAGS})


# Option.
if(CMAKE_BUILD_TYPE STREQUAL Debug)
    target_compile_options(MsaCalc PRIVATE -DDEBUG -g -O0 -Wall -fopenmp)
    set(CMAKE_DEBUG_POSTFIX "_d")
endif(CMAKE_BUILD_TYPE STREQUAL Debug)

if(CMAKE_BUILD_TYPE STREQUAL Release)
    target_compile_options(MsaCalc PRIVATE -g -O3 -Wall -fopenmp -march=native)
endif(CMAKE_BUILD_TYPE STREQUAL Release)


file(COPY ${CMAKE_SOURCE_DIR}/Test DESTINATION .)
