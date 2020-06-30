#include "Simulation.h"

#include <iostream>
#include <cassert>
#include <insitu/sensei/multiGrid.h>
#include <insitu/sensei/unstructuredGrid.h>
#include <core/unstr.h>

using namespace vistle::insitu::sensei;
using namespace vistle::insitu;
using std::cerr; using std::endl;
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
vistle::insitu::sensei::Grid Simulation::getMesh(const std::string& name) const
{
	cerr << "preparing mesh " << name << endl;
	if (name != meshName)
	{
		cerr << "wrong mesh name " << name << endl;
		return nullptr;
	}




	auto multiMesh = std::make_unique<MultiGrid>( name );
	std::array<double, 2> xOffsets{ 0,1 }, yOffsets{ 1,0 };
	for (size_t i = 0; i < 2; i++)
	{
		auto subMultiMesh = std::make_unique<MultiGrid>(name);
		subMultiMesh->addGrid(createMesh(name, xOffsets[i] * 10, yOffsets[i] *10, 0));
		subMultiMesh->addGrid(createMesh(name, xOffsets[i] * 5, yOffsets[i] * 5, 1));
		multiMesh->addGrid(std::move(subMultiMesh));

	}


	return std::move(multiMesh);

}

Grid Simulation::createMesh(const std::string& name, double xOffset, double yOffset, size_t index) const
{
	std::array<double, 3> dist{ 1,2,1 };
	std::array<Array, 3> arrays;

	for (size_t i = 0; i < 3; i++)
	{
		arrays[i] = Array::create<double>(name, dims[0] * dims[1] * dims[2], DataMapping::Vertex);
	}
	double* x = arrays[0].dataAs<double>();
	double* y = arrays[1].dataAs<double>();
	double* z = arrays[2].dataAs<double>();

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
	auto grid = std::make_unique<vistle::insitu::sensei::UnstructuredMesh>(name, std::move(arrays), dims);

	size_t numElem = 1;
	for (size_t i = 0; i < 3; i++)
	{
		numElem *= std::max(size_t{ 1 }, dims[i]);
	}
	size_t numCellVert = 8;
	auto clArray = Array::create<vistle::Index>(name, numElem * numCellVert);
	auto elArray = Array::create<vistle::Index>(name, numElem);
	auto tlArray = Array::create<vistle::Byte>(name, numElem);

	vistle::Index elem = 0;
	vistle::Index idx = 0;
	vistle::Index* cl = clArray.dataAs<vistle::Index>();
	vistle::Index* el = elArray.dataAs<vistle::Index>();
	vistle::Byte* tl = tlArray.dataAs<vistle::Byte>();

	/*unsigned type = vistle::UnstructuredGrid::POINT;
	if (dim[2] > 1) {
		type = UnstructuredGrid::HEXAHEDRON;
	}
	else if (dim[1] > 1) {
		type = UnstructuredGrid::QUAD;
	}
	else if (dim[0] > 1) {
		type = UnstructuredGrid::BAR;
	}

	for (Index ix = 0; ix < nx; ++ix) {
		for (Index iy = 0; iy < ny; ++iy) {
			for (Index iz = 0; iz < nz; ++iz) {
				cl[idx++] = UniformGrid::vertexIndex(ix, iy, iz, dim);
				if (dim[0] > 1) {
					cl[idx++] = UniformGrid::vertexIndex(ix + 1, iy, iz, dim);
					if (dim[1] > 1) {
						cl[idx++] = UniformGrid::vertexIndex(ix + 1, iy + 1, iz, dim);
						cl[idx++] = UniformGrid::vertexIndex(ix, iy + 1, iz, dim);
						if (dim[2] > 1) {
							cl[idx++] = UniformGrid::vertexIndex(ix, iy, iz + 1, dim);
							cl[idx++] = UniformGrid::vertexIndex(ix + 1, iy, iz + 1, dim);
							cl[idx++] = UniformGrid::vertexIndex(ix + 1, iy + 1, iz + 1, dim);
							cl[idx++] = UniformGrid::vertexIndex(ix, iy + 1, iz + 1, dim);
						}
					}
				}

				tl[elem] = type;
				if ((ix < ghostWidth[0][0] || ix + ghostWidth[0][1] >= nx)
					|| (iy < ghostWidth[1][0] || iy + ghostWidth[1][1] >= ny)
					|| (iz < ghostWidth[2][0] || iz + ghostWidth[2][1] >= nz)) {
					tl[elem] |= UnstructuredGrid::GHOST_BIT;
				}

				++elem;
				el[elem] = idx;
			}
		}
	}*/




	return grid;
}

vistle::insitu::Array Simulation::getVar(const std::string& name) const
{
	cerr << "preparing variable data " << name << endl;
	if (name != varname)
	{
		cerr << "wrong variable name " << name << endl;
		return Array{};
	}
	auto multi = Array::create<Array>(name, 2, DataMapping::Vertex);
	for (size_t j = 0; j < 2; j++)
	{
		auto subArray = Array::create<Array>(name, 2, DataMapping::Vertex);
		for (size_t i = 0; i < 2; i++)
		{

			subArray.dataAs<Array>()[i] = Array::create<float>(name, dims[0] * dims[1] * dims[2], DataMapping::Vertex);
			for (size_t j = 0; j < subArray.dataAs<Array>()[i].size(); j++)
			{
				float* data = subArray.dataAs<Array>()[i].dataAs<float>();
				data[j] = j * timestep + 3 * i;
			}
		}
		multi.dataAs<Array>()[j] = std::move(subArray);

	}
	return multi;
}

size_t Simulation::nextStep()
{
	return timestep++;
}
