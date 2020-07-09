#include "Simulation.h"

#include <iostream>
#include <cassert>
#include <insitu/sensei/multiGrid.h>
#include <insitu/sensei/unstructuredGrid.h>
#include <insitu/core/transformArray.h>

#include <core/unstr.h>
#include <core/archives.h>
#include <core/uniformgrid.h>
#include <core/structuredgrid.h>
#include <core/rectilineargrid.h>

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


vistle::Object::ptr Simulation::makeUnstr(const std::string& name)
{
	using namespace vistle;
	const size_t nx = std::max(dims[0] - 1, size_t(1));
	const size_t ny = std::max(dims[1] - 1, size_t(1));
	const size_t nz = std::max(dims[2] - 1, size_t(1));
	size_t numElem = nx * ny * nz;
	size_t numCellVert = 8;
	auto clArray = Array::create<vistle::Index>(name, numElem * numCellVert);
	auto elArray = Array::create<vistle::Index>(name, numElem);
	auto tlArray = Array::create<vistle::Byte>(name, numElem);

	vistle::Index elem = 0;
	vistle::Index idx = 0;
	vistle::Index* cl = clArray.dataAs<vistle::Index>();
	vistle::Index* el = elArray.dataAs<vistle::Index>();
	vistle::Byte* tl = tlArray.dataAs<vistle::Byte>();

	unsigned type = UnstructuredGrid::HEXAHEDRON;
	return nullptr;

	//for (size_t ix = 0; ix < nx; ++ix) {
	//	for (size_t iy = 0; iy < ny; ++iy) {
	//		for (size_t iz = 0; iz < nz; ++iz) {
	//			cl[idx++] = UniformGrid::vertexIndex(ix, iy, iz, dims.data());
	//			if (dims[0] > 1) {
	//				cl[idx++] = UniformGrid::vertexIndex(ix + 1, iy, iz, dims.data());
	//				if (dims[1] > 1) {
	//					cl[idx++] = UniformGrid::vertexIndex(ix + 1, iy + 1, iz, dims.data());
	//					cl[idx++] = UniformGrid::vertexIndex(ix, iy + 1, iz, dims.data());
	//					if (dims[2] > 1) {
	//						cl[idx++] = UniformGrid::vertexIndex(ix, iy, iz + 1, dims.data());
	//						cl[idx++] = UniformGrid::vertexIndex(ix + 1, iy, iz + 1, dims.data());
	//						cl[idx++] = UniformGrid::vertexIndex(ix + 1, iy + 1, iz + 1, dims.data());
	//						cl[idx++] = UniformGrid::vertexIndex(ix, iy + 1, iz + 1, dims.data());
	//					}
	//				}
	//			}

	//			tl[elem] = type;
	//			/*				if ((ix < ghostWidth[0][0] || ix + ghostWidth[0][1] >= nx)
	//								|| (iy < ghostWidth[1][0] || iy + ghostWidth[1][1] >= ny)
	//								|| (iz < ghostWidth[2][0] || iz + ghostWidth[2][1] >= nz)) {
	//								tl[elem] |= UnstructuredGrid::GHOST_BIT;
	//							*/
	//		}

	//		++elem;
	//		el[elem] = idx;
	//	}
	//}

	//auto grid = std::make_unique<insitu::sensei::UnstructuredMesh>(name);
	//grid->cl = std::move(clArray);
	//grid->el = std::move(elArray);
	//grid->tl = std::move(tlArray);
	//grid->clToVistle = [](Index* vistleCL, const Array& myCL)->bool {
	//	vistle::insitu::transformArray(myCL, vistleCL);
	//	return true; };
	//grid->elToVistle = [](Index* vistleCL, const Array& myCL)->bool {
	//	vistle::insitu::transformArray(myCL, vistleCL);
	//	return true; };
	//grid->tlToVistle = [](Byte* vistleCL, const Array& myCL)->bool {
	//	vistle::insitu::transformArray(myCL, vistleCL);
	//	return true; };
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
