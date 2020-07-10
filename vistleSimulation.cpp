#include "Simulation.h"

#include <mpi.h>
#include <iostream>


using namespace vistle::insitu;

//struct Test {
//	Test():t(new int) {
//		*t.get() = 5;
//	};
//	std::shared_ptr<int> t;
//	std::shared_ptr<int> add() {
//		(*t.get())++;
//		return t;
//	}
//};
//
//struct Run {
//	Run() {
//		std::function<std::shared_ptr<int>()> f;
//		Test t;
//
//		std::cerr << "t.value = " << *t.t.get() << std::endl;
//
//		f = std::bind(&Test::add, t);
//		f();
//		std::cerr << "t.value = " << *t.t.get() << std::endl;
//	}
//
//
//};

int main(int argc, char* argv[]){

	MPI_Init(&argc, &argv);
	
	Simulation sim;
	sim.run();
	MPI_Finalize();

}