#pragma once
#include <insitu/sensei/senseiInterface.h>
#include <insitu/sensei/structuredGrid.h>
class Simulation
{
public:
	void intitMesh();
	vistle::insitu::sensei::Grid getMesh(const std::string& name) const;
	vistle::insitu::sensei::Grid createMesh(const std::string& name, double xOffset, double yOffset,  size_t index) const;
	vistle::insitu::Array getVar(const std::string& name) const;

	size_t nextStep();
	std::string meshName = "mesh", varname = "data";
private:
	constexpr static std::array<size_t, 3> dims{ 2,4,6 };
	size_t timestep = 0;
	std::array<float, dims[0]> xCoords;
	std::array<float, dims[1]> yCoords;
	std::array<float, dims[2]> zCoords;
	std::array<float*, 3> coords{ xCoords.data(), yCoords.data(), zCoords.data() };
};

