/*
* @file 	waterentry_elastic_case.h
*/

#ifndef	TANK_CASE_H
#define TANK_CASE_H

#include "sphinxsys.h"
#define PI (3.14159265358979323846)
using namespace SPH;

/*@brief Basic geometry parameters and numerical setup.
*/

Real resolution_ref = 0.0075;   /* Initial particle spacing*/


/* Domain bounds of the system*/
BoundingBox system_domain_bounds(Vec3d(-0.2, -0.05,-0.2), Vec3d(0.2, 1.0,0.2));


/*
Material properties of the fluid.
*/
Real rho0_f = 1000.0;         /*Fluid density*/
Real rho0_a = 1.226;            /*Air density*/
Real rho0_s = 7890;
Real gravity_g = 9.81;        /*Gravity force of fluid*/
Real U_f = 2.0* sqrt(gravity_g * 0.5);	/**< Characteristic velocity. */
Real U_g = 2.0* sqrt(gravity_g * 0.5);  	/**< dispersion velocity in shallow water. */
Real c_f = 10.0 * SMAX(U_g, U_f);	/**< Reference sound speed. */
Real f = 1.0;
Real a = 0.08;
Real poisson = 0.27; 		/**< Poisson ratio.*/
Real Ae = 135.0e9; 			/**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae;
Real mu_water = 653.9e-6;
Real mu_air = 20.88e-6;
Real length_scale = 1.0;
Vec3d translation(0, 0.0, 0);

std::string fuel_tank_outer = "./input/tank_outer.STL";
std::string fuel_tank_inner = "./input/tank_inner.STL";
std::string water_05 = "./input/water_05.STL";
std::string air_05 = "./input/gas_05.STL";
std::string probe_s1_shape = "./input/ProbeS1.STL";
std::string probe_s2_shape = "./input/ProbeS2.STL";
std::string probe_s3_shape = "./input/ProbeS3.STL";

/*
Fuel Tank.
*/
class Tank : public ComplexShape
{
public:
	explicit Tank(const std::string &shape_name) :ComplexShape(shape_name)
	{
		/** Geometry definition. */

		add<TriangleMeshShapeSTL>(fuel_tank_outer, translation, length_scale,"OuterWall");
		subtract<TriangleMeshShapeSTL>(fuel_tank_inner, translation, length_scale,"InnerWall");
	}
};
class WaterBlock : public ComplexShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(water_05, translation, length_scale);
	}
};

class AirBlock : public ComplexShape
{
public:
	explicit AirBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(air_05, translation, length_scale);
	}
};

class VariableGravity : public Gravity
{
	Real time_ = 0;
public:
	VariableGravity() : Gravity(Vecd(0.0, -gravity_g, 0.0)) {};
	virtual Vecd InducedAcceleration(Vecd& position) override
	{
		time_= GlobalStaticVariables::physical_time_;
		global_acceleration_[0] = 4.0 * PI * PI * f * f * a * sin(2 * PI * f * time_);
		return global_acceleration_;
	}
};
#endif //TANK_CASE_H