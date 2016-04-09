cmake_minimum_required (VERSION 2.6)
project (findSpg)

add_subdirectory(spglib)

set(FIND_SPG_SRCS
    src/findSpg.cpp
    src/crystal.cpp
    src/elemInfo.cpp)

include_directories(${findSpg_SOURCE_DIR}/include)

if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
endif(UNIX)

add_executable (findSpg ${FIND_SPG_SRCS})
target_link_libraries (findSpg spglib)
