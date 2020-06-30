cmake_minimum_required(VERSION 3.14)
project(VistleSImulation)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)




find_package(vistle_sensei REQUIRED)
find_package(MPI REQUIRED)
add_executable(VistleSimulation vistleSimulation.cpp Simulation.h Simulation.cpp)

target_link_libraries(VistleSimulation 
	PRIVATE Vistle::vistle_sensei
	PRIVATE MPI::MPI_C
	)

message("BOOST_LOCATION = " ${BOOST_LOCATION})
message("BOOST_INCLUDE_DIRECTORIES = " ${BOOST_INCLUDE_LOCATION})