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

Real resolution_ref = 0.005;   /* Initial particle spacing*/


/* Domain bounds of the system*/
BoundingBox system_domain_bounds(Vec3d(-0.3, -0.3,-0.3), Vec3d(0.3, 0.5,0.3));


/*
Material properties of the fluid.
*/
Real rho0_f = 1000.0;         /*Fluid density*/
Real rho0_a = 1.0;            /*Air density*/
Real gravity_g = 9.81;        /*Gravity force of fluid*/
Real U_f = 2.0* sqrt(gravity_g * 0.174);	/**< Characteristic velocity. */
Real U_g = 2.0* sqrt(gravity_g * 0.174);  	/**< dispersion velocity in shallow water. */
Real c_f = 10.0 * SMAX(U_g, U_f);	/**< Reference sound speed. */


Real length_scale = 1.0;
Vec3d translation(0, 0.175, 0);

std::string fuel_tank_outer = "./input/validation_tank_outer_slim.STL";
std::string fuel_tank_inner = "./input/validation_tank_inner.STL";
std::string water_05 = "./input/validation_water.STL";
std::string air_05 = "./input/validation_air.STL";
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
		if (time_ > 0.25)
		{
			global_acceleration_[0] = 4 * PI*PI* 1.63*1.63*0.0075*sin(2 * PI*1.63*(time_-0.25));
		}
		
		return global_acceleration_;
	}
};

class ProbeS1 : public ComplexShape
{
public:
	explicit ProbeS1(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vec3d translation_probe(0.0, 0.0, 0.0);
		add<TriangleMeshShapeSTL>(probe_s1_shape, translation_probe, length_scale);
	}
};

class ProbeS2 : public ComplexShape
{
public:
	explicit ProbeS2(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vec3d translation_probe_2(0.0, 0.0, 0.0);
		add<TriangleMeshShapeSTL>(probe_s2_shape, translation_probe_2, length_scale);
	}
};

class ProbeS3 : public ComplexShape
{
public:
	explicit ProbeS3(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vec3d translation_probe_3(0.0,0.0, 0.0);
		add<TriangleMeshShapeSTL>(probe_s3_shape, translation_probe_3, length_scale);
	}
};

/*
class SloshMaking : public solid_dynamics::ConstrainSolidBodyRegion
{

	Real wave_freq_;
	Real wave_stroke_;
	Real time_;

	virtual Vecd getDisplacement(Vecd &pos_0, Vecd &pos_n) override
	{
		Vecd displacement(0);
		displacement[0] = -(wave_stroke_ / (2 * PI * wave_freq_ * 2 * PI * wave_freq_)) * sin(2 * PI * wave_freq_ * time_);
		return pos_0 + displacement;
	}

	virtual Vecd getVelocity(Vecd &pos_0, Vecd &pos_n, Vecd &vel_n) override
	{
		Vecd velocity(0);
		velocity[0] = - (wave_stroke_ /( 2*PI*wave_freq_)) * cos(2 * PI*wave_freq_ * time_);
		return velocity;
	}

	virtual Vecd getAcceleration(Vecd &pos_0, Vecd &pos_n, Vecd &dvel_dt) override
	{
		Vecd acceleration(0);
		acceleration[0] =  wave_stroke_ *  sin(2*PI*wave_freq_ * time_);
		return acceleration;
	}

	virtual void setupDynamics(Real dt = 0.0) override
	{
		body_->setNewlyUpdated();
		time_ = GlobalStaticVariables::physical_time_;
	}

	void computeWaveStrokeAndFrequency()
	{

		wave_stroke_ = 0.5;
		if (time_ < 9.5 || (time_ > 18.0 && time_ <= 28.0))
		{
			wave_freq_=1.4;
		}

		if (time_ > 28.0)
		{
			wave_freq_ = -0.0000000597 *pow(time_, 5) + 0.0000046039*pow(time_, 4)
				+ 0.0000486365*pow(time_, 3) - 0.0129251176 *pow(time_, 2)
				+ 0.3917360229*time_ - 2.3043720050;
		}

		if (time_ >= 9.5 && time_ <= 18.0)
		{
			wave_freq_ = 0.0000001901*pow(time_, 5) - 0.0000264917*pow(time_, 4)
				+ 0.001097686*pow(time_, 3) - 0.0196394691*pow(time_, 2)
				+ 0.1589433209*time_ + 0.9219187406;
		}

		std::cout << "Wave stroke: " << wave_stroke_ << " Wave frequency: " << wave_freq_ << std::endl;
	}

public:
	SloshMaking(SolidBody &solid_body, BodyPartByParticle &constrained_region)
		: ConstrainSolidBodyRegion(solid_body, constrained_region), time_(0.0)
	{
		computeWaveStrokeAndFrequency();
	}
};
*/
#endif //TANK_CASE_H