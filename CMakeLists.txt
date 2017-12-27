#to collect compiler flags and specific system settings to compile the TopologyTracer2D3D
cmake_minimum_required(VERSION 2.6)
project(topotracer2d3d)

##define level of code optimization, -O0 nothing (for debug), -O3 aggressive for production jobs
set(CMAKE_BUILD_DIR "build")
set(OPTLEVEL "-O0")
##choose compiler
set(TOPOTRACER_EMPLOYING_INTELCOMPILER ON)
set(TOPOTRACER_EMPLOYING_GNUCOMPILER OFF)


##utilize TetGen for 3d grain hull tetrahedralization
##this feature is BETA, consult source code before switching on!
set(TOPOTRACER_EMPLOYING_TETGEN OFF)

##end of user-interaction
#################################################################################################################
message([STATUS] "We utilize optlevel ${OPTLEVEL}")

if(TOPOTRACER_EMPLOYING_INTELCOMPILER AND TOPOTRACER_EMPLOYING_GNUCOMPILER)
	message([FATAL_ERROR] "You cannot utilize two compiler at the same time!")
endif()

if(TOPOTRACER_EMPLOYING_INTELCOMPILER)
	message([STATUS] "Employing the Intel compiler!")
	set(MPILARGEFILES "-D_FILE_OFFSET_BITS=64")
	set(OOMP "${OPTLEVEL} -qopenmp -lpthread")
	add_definitions("-std=c++0x")
	add_definitions("-Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wpointer-sign -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-variable") #-Wunused-function -Wunused-but-set-variable -Wunused-variable")
elseif(TOPOTRACER_EMPLOYING_GNUCOMPILER)
	message([STATUS] "Employing the GNU compiler!")
	set(OOMP "${OPTLEVEL} -fopenmp -lpthread")
	add_definitions("-std=c++11")
	add_definitions("-Wall -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable")
else()
	message([FATAL_ERROR] "You have to utilize a compiler!")
endif()

##Compiler flag passing
set(CMAKE_C_FLAGS "${OPTLEVEL} ${OOMP}") 
set(CMAKE_CXX_FLAGS "${OPTLEVEL} ${OOMP}") 

##assure MPI is there to employ process parallelism
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

##add source code snippets
add_executable(topotracer2d3d
			../src/thirdparty/PolyTri/predicates.cc 
			../src/thirdparty/PolyTri/PolyTri_Geometry.cpp
			../src/TopologyTracer2D3D_Allocator.cpp
			../src/TopologyTracer2D3D_Settings.cpp
			../src/TopologyTracer2D3D_Math.cpp
			../src/TopologyTracer2D3D_Topohdl.cpp
			../src/TopologyTracer2D3D_Main.cpp
)

##assure the Boost package is there to allow and ease the reading of files and folders from the filesystem
message([STATUS] "We need Boost as well...")
find_package(Boost COMPONENTS system filesystem REQUIRED)

##additional TopoTracer functionality
if(TOPOTRACER_EMPLOYING_TETGEN)
	message([STATUS] "TetGen is activated!")
	set(TETGEN_FLAGS "-L/home/mk958655/2016TopologyTracer2D3D/src/ -ltet")
else()
	message([STATUS] "TetGen is not activated!")
	set(TETGEN_FLAGS "")
endif()

##linking the program together
target_link_libraries(topotracer2d3d ${OOMP} ${TETGEN_FLAGS} -lm -lnuma ${MPI_LIBRARIES} ${MPILARGEFILES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} )

##MPI compiler flags
if(MPI_COMPILE_FLAGS)
  set_target_properties(topotracer2d3d PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(topotracer2d3d PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
