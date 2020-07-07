#include "Simulation.h"

#include <insitu/sensei/senseiInterface.h>
#include <insitu/sensei/gridInterface.h>

#include <mpi.h>
#include <iostream>


using namespace vistle::insitu;

struct Test {
	int t = 0;
	void add() {
		t++;
	}
};

int main(int argc, char* argv[]){

	Test t;
	std::function<void()>add = std::bind(&Test::add, &t);
	for (size_t i = 0; i < 3; i++)
	{
		add();
	}
	std::cerr << " t.t = " << t.t << std::endl;
	MPI_Init(&argc, &argv);
	std::cerr << "hello world" << std::endl;
	Simulation sim;
	sim.intitMesh();
	vistle::insitu::MetaData meta{};
	auto it = meta.addMesh(sim.meshName);
	meta.addVariable(sim.varname, it);
	std::function<vistle::insitu::sensei::Grid(const std::string&)> m = std::bind(&Simulation::getMesh, &sim, std::placeholders::_1);
	std::function<vistle::insitu::Array(const std::string&)> v = std::bind(&Simulation::getVar, &sim, std::placeholders::_1);
	auto vistleSensei = vistle::insitu::sensei::createSenseiInterface(true, 0, 1, std::move(meta), vistle::insitu::sensei::Callbacks{ m, v });
	std::vector<Array> vars{};
	while (vistleSensei->Execute(sim.nextStep()))
	{

	}
	vistleSensei->Finalize();
	MPI_Finalize();

}