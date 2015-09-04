#to compile a program for tracing simulation step resolved single-grain evolution
cmake_minimum_required(VERSION 2.6)
project(topologytracer)

set ( CMAKE_CXX_FLAGS -O2 )
set ( CMAKE_C_FLAGS -O2 )



#larger than 2^32 bit files on 64bit architecture?
#compile with option -D_FILE_OFFSET_BITS=64

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
add_executable(topologytracer
			../src/topologytracer_math.cpp
			../src/topologytracer_random.cpp
			
			../src/topologytracer.cpp
			../src/topologytracer_kernel.cpp
)
target_link_libraries(topologytracer ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(topologytracer PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(topologytracer PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()