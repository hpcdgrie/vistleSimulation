#pragma once

#include <insitu/sensei/sensei.h>

class Simulation
{
public:
	Simulation();
	void run();
	~Simulation();
	void intitMesh();
	vistle::Object::ptr getMesh(const std::string& name) const;
	vistle::Object::ptr createMesh(const std::string& name, double xOffset, double yOffset,  size_t index) const;
	vistle::Object::ptr makeUnstr(const std::string& name);
	vistle::DataBase::ptr getVar(const std::string& name) const;

	size_t nextStep();
	std::string meshName = "mesh", varname = "data";
private:
	constexpr static std::array<size_t, 3> dims{ 2,4,6 };
	size_t timestep = 0;
	std::array<vistle::Index, dims[0]> xCoords;
	std::array<vistle::Index, dims[1]> yCoords;
	std::array<vistle::Index, dims[2]> zCoords;
	std::array<vistle::Index*, 3> coords{ xCoords.data(), yCoords.data(), zCoords.data() };
	std::unique_ptr<vistle::insitu::sensei::SenseiAdapter> m_adapter;
};

