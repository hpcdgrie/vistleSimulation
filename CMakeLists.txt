cmake_minimum_required(VERSION 3.14)
project(VistleSImulation)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)


add_definitions(-DBOOST_MPI_SOURCE)


find_package(vistle_sensei REQUIRED)
find_package(MPI REQUIRED)
add_executable(VistleSimulation vistleSimulation.cpp Simulation.h Simulation.cpp)

target_link_libraries(VistleSimulation 
	PRIVATE Vistle::vistle_sensei
	PRIVATE MPI::MPI_C
	)
include_directories("H:/src/vistle/3rdparty/yas/include")
message("BOOST_LOCATION = " ${BOOST_LOCATION})
message("BOOST_INCLUDE_DIRECTORIES = " ${BOOST_INCLUDE_LOCATION})