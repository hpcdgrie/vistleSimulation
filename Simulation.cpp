#include "Simulation.h"

#include <core/unstr.h>
#include <core/uniformgrid.h>
#include <core/structuredgrid.h>
#include <core/rectilineargrid.h>

#include <iostream>
#include <cassert>
#include <algorithm>

using namespace vistle::insitu::sensei;
using namespace vistle::insitu;
using std::cerr; using std::endl;
Simulation::Simulation()
{
	intitMesh();
	vistle::insitu::MetaData meta{};
	auto it = meta.addMesh(meshName);
	meta.addVariable(varname, it);


	//std::function<vistle::Object::ptr(const std::string&)> gm = [=](const std::string& name) {
	//	return this->getMesh(name);
	//};
	//std::function<vistle::DataBase::ptr(const std::string&)> gv = [=](const std::string& name) {
	//	return this->getVar(name);
	//};
	std::function<vistle::Object::ptr(const std::string&)> m = std::bind(&Simulation::getMesh, this, std::placeholders::_1);
	std::function<vistle::DataBase::ptr(const std::string&)> v = std::bind(&Simulation::getVar, this, std::placeholders::_1);
	m_adapter = std::make_unique<vistle::insitu::sensei::SenseiAdapter>(true, 0, 1, std::move(meta), vistle::insitu::sensei::Callbacks{  m, v });


}

void Simulation::run() {

	cerr << "start running" << endl;
	while (m_adapter->Execute(nextStep()))
	{

	}
}

Simulation::~Simulation() {
	m_adapter->Finalize();
}
void Simulation::intitMesh()
{
	for (size_t i = 0; i < 3; i++)
	{
		for (size_t j = 0; j < dims[i]; j++)
		{
			coords[i][j] = j;
		}
	}
}
template<typename T>
T vertexIndex(const T ix, const T iy, const T iz, const std::array<T, 3> dims) {
	assert(ix < dims[0]);
	assert(iy < dims[1]);
	assert(iz < dims[2]);
	return (ix * dims[1] + iy) * dims[2] + iz;
}
vistle::Object::ptr Simulation::getMesh(const std::string& name) const
{
	cerr << "preparing mesh " << name << endl;
	if (name != meshName)
	{
		cerr << "wrong mesh name " << name << endl;
		return nullptr;
	}


	return createMesh(name, 0, 0, 0);

}

vistle::Object::ptr Simulation::createMesh(const std::string& name, double xOffset, double yOffset, size_t index) const
{
	//return nullptr;
	using namespace vistle;
	std::array<double, 3> dist{ 1,2,1 };


	auto grid = m_adapter->createVistleObject<vistle::StructuredGrid>(dims[0], dims[1], dims[2]);
	vistle::Scalar* x = grid->x().data();
	vistle::Scalar* y = grid->y().data();
	vistle::Scalar* z = grid->z().data();

	//#pragma omp parallel for
	for (size_t i = 0; i < dims[0]; ++i) {
		//#pragma omp parallel for
		for (size_t j = 0; j < dims[1]; ++j) {
			//#pragma omp parallel for
			for (size_t k = 0; k < dims[2]; ++k) {
				auto idx = vertexIndex(i, j, k, dims);
				assert(idx < dims[0] * dims[1] * dims[2]);
				x[idx] = xOffset + i * dist[0];
				y[idx] = yOffset + j * dist[1];
				z[idx] = k * dist[2];
			}
		}
	}
	return grid;
}


vistle::DataBase::ptr Simulation::getVar(const std::string& name) const
{
	cerr << "preparing variable data " << name << endl;
	if (name != varname)
	{
		cerr << "wrong variable name " << name << endl;
		return nullptr;
	}
	auto var = m_adapter->createVistleObject<vistle::Vec<vistle::Scalar>>(dims[0] * dims[1] * dims[2]);
	for (size_t j = 0; j < var->getSize(); j++)
	{
		var->x().data()[j] = j;
	}
	return var;
}

size_t Simulation::nextStep()
{
	return timestep++;
}
