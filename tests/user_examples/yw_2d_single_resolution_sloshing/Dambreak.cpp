/**
 * @file 	dambreak.cpp
 * @brief 	2D dambreak example.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for fluid simulation.
 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
#define PI (3.14159265358979323846)
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.08;					/**< Water tank length. */
Real DH = 0.12;					/**< Water tank height. */
Real LL = 0.08;						/**< Water column length. */
Real LH = 0.06;						/**< Water column height. */

Real particle_spacing_ref = 0.002;	/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Thickness of tank wall. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 818.0;					 /**< Reference density of fluid. */
Real rho0_v = 0.899;						 /**< Reference density of air. */
Real gravity_g = 9.81;					 /**< Gravity. */
Real U_max = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;				 /**< Reference sound speed. */
Real diffusion_coff_liquid = 6.1125e-8;
Real diffusion_coff_vapor = 26.7468e-6;
Real mu_f = 818e-6;
Real mu_v = 1.7894e-5;

//----------------------------------------------------------------------
//	Geometric elements used in shape modeling.
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(0.0, 0.0));
	water_block_shape.push_back(Vecd(0.0, LH));
	water_block_shape.push_back(Vecd(LL, LH));
	water_block_shape.push_back(Vecd(LL, 0.0));
	water_block_shape.push_back(Vecd(0.0, 0.0));
	return water_block_shape;
}

std::vector<Vecd> createOuterWallShape()
{
	std::vector<Vecd> outer_wall_shape;
	outer_wall_shape.push_back(Vecd(-BW, -BW));
	outer_wall_shape.push_back(Vecd(-BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, -BW));
	outer_wall_shape.push_back(Vecd(-BW, -BW));

	return outer_wall_shape;
}

std::vector<Vecd> createInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(0.0, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, DH));
	inner_wall_shape.push_back(Vecd(DL, DH));
	inner_wall_shape.push_back(Vecd(DL, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, 0.0));

	return inner_wall_shape;
}
//----------------------------------------------------------------------
//	cases-dependent geometric shape for water block.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
	}
};

//----------------------------------------------------------------------
//	cases-dependent geometric shape for air block.
//----------------------------------------------------------------------
class VaporBlock : public MultiPolygonShape
{
public:
	explicit VaporBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::sub);
	}
};
//----------------------------------------------------------------------
//	Wall boundary shape definition.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
public:
	explicit WallBoundary(const std::string& shape_name) : ComplexShape(shape_name)
	{
		add<MultiPolygonShape>(MultiPolygon(createOuterWallShape()), "OuterWall");
		subtract<MultiPolygonShape>(MultiPolygon(createInnerWallShape()), "InnerWall");
	}
};


class VariableGravity : public Gravity
{
	Real time_ = 0;
public:
	VariableGravity() : Gravity(Vecd(0.0, -gravity_g)) {};
	virtual Vecd InducedAcceleration(Vecd& position) override
	{
		time_ = GlobalStaticVariables::physical_time_;
		//if (time_ <= 0.24)
		//{
			//global_acceleration_[0] = 1312.547554*pow(time_, 4) - 449.087146*pow(time_, 3) - 79.178542*pow(time_, 2) - 6.428713*time_ + 7.956987;
		global_acceleration_[0] = 0.6 * gravity_g * cos(2 * PI * time_);
		//}
		//else
		//{
			//global_acceleration_[0] = 0.0;
		//}
		return global_acceleration_;
	}
};

Real densed_point = 0.06;


Real h = 1.3 * particle_spacing_ref;
MultiPolygon createWaveProbeShape4()
{
	std::vector<Vecd> pnts;
	pnts.push_back(Vecd(0.4 - h, 0.0));
	pnts.push_back(Vecd(0.4 - h, 0.6));
	pnts.push_back(Vecd(0.4 + h, 0.6));
	pnts.push_back(Vecd(0.4 + h, 0.0));
	pnts.push_back(Vecd(0.4 - h, 0.0));

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
	return multi_polygon;
}




//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char* av[])
{
	//----------------------------------------------------------------------
	//	Build up an SPHSystem.
	//----------------------------------------------------------------------
	BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
	GlobalStaticVariables::physical_time_ = 0.0;
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	sph_system.run_particle_relaxation_ = false;
	/** Tag for computation start with relaxed body fitted particles distribution. */
	sph_system.reload_particles_ = false;
	/* Tag for computation from restart files. 0: start with initial condition. */
	sph_system.restart_step_ = 0;
#ifdef BOOST_AVAILABLE
	sph_system.handleCommandlineOptions(ac, av);
#endif

	IOEnvironment io_environment(sph_system);


	//----------------------------------------------------------------------
	//	Creating bodies with corresponding materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
	//water_block.defineAdaptation<ParticleRefinementByPosition>(1.3, 1.0, 2, densed_point);
	water_block.defineAdaptation<SPHAdaptation>(1.3, 2);
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
	//water_block.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	(!sph_system.run_particle_relaxation_ && sph_system.reload_particles_)
		? water_block.generateParticles<ParticleGeneratorReload>(io_environment, water_block.getName())
		: water_block.generateParticles<ParticleGeneratorLattice>();
	//water_block.addBodyStateForRecording<Real>("SmoothingLengthRatio");

	/*FluidBody vapor_block(sph_system, makeShared<VaporBlock>("VaporBody"));
	vapor_block.defineAdaptation<ParticleRefinementByPosition>(1.3, 1.0, 1, densed_point);
	vapor_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_v, c_f);
	vapor_block.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	(!sph_system.run_particle_relaxation_ && sph_system.reload_particles_)
		? vapor_block.generateParticles<ParticleGeneratorReload>(io_environment, vapor_block.getName())
		: vapor_block.generateParticles<ParticleGeneratorMultiResolutionByPosition>();
*/
	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
	//wall_boundary.defineAdaptation<SPHAdaptation>(1.15, 1.0);
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	//wall_boundary.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	/*(!sph_system.run_particle_relaxation_ && sph_system.reload_particles_)
		? wall_boundary.generateParticles<ParticleGeneratorReload>(io_environment, wall_boundary.getName())
		: */wall_boundary.generateParticles<ParticleGeneratorLattice>();
		wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");


		//----------------------------------------------------------------------
		//	Define body relation map.
		//	The contact map gives the topological connections between the bodies.
		//	Basically the the range of bodies to build neighbor particle lists.
		//----------------------------------------------------------------------
		//----------------------------------------------------------------------
		//	Define body relation map.
		//	The contact map gives the topological connections between the bodies.
		//	Basically the the range of bodies to build neighbor particle lists.
		//----------------------------------------------------------------------

		//InnerRelation wall_boundary_inner(wall_boundary);
		InnerRelation water_inner(water_block);
		//AdaptiveInnerRelation vapor_inner(vapor_block);
		if (sph_system.run_particle_relaxation_)
		{
			//----------------------------------------------------------------------
			//	Methods used for particle relaxation.
			//----------------------------------------------------------------------
			/** Random reset the insert body particle position. */
			//SimpleDynamics<RandomizeParticlePosition> random_wall_boundary_particles(wall_boundary);
			SimpleDynamics<RandomizeParticlePosition> random_water_particles(water_block);
			//SimpleDynamics<RandomizeParticlePosition> random_vapor_particles(vapor_block);
			/** Write the body state to Vtp file. */
			//BodyStatesRecordingToVtp write_wall_boundary_to_vtp(io_environment, { &wall_boundary });
			BodyStatesRecordingToVtp write_water_to_vtp(io_environment, { &water_block });
			//BodyStatesRecordingToVtp write_vapor_to_vtp(io_environment, { &vapor_block });
			/** Write the particle reload files. */
			//ReloadParticleIO write_wall_boundary_particle_reload_files(io_environment, wall_boundary, "Wall_Boundary");
			ReloadParticleIO write_water_particle_reload_files(io_environment, water_block, "WaterBody");
			//ReloadParticleIO write_vapor_particle_reload_files(io_environment, vapor_block, "VaporBody");
			/** A  Physics relaxation step. */
			relax_dynamics::RelaxationStepInner water_relaxation_step_inner(water_inner);
			//relax_dynamics::RelaxationStepInner vapor_relaxation_step_inner(vapor_inner);
			//relax_dynamics::RelaxationStepInner wall_boundary_relaxation_step_inner(wall_boundary_inner);
			SimpleDynamics<relax_dynamics::UpdateSmoothingLengthRatioByPosition> update_smoothing_length_ratio_water(water_block);
			//SimpleDynamics<relax_dynamics::UpdateSmoothingLengthRatioByPosition> update_smoothing_length_ratio_vapor(vapor_block);
			//----------------------------------------------------------------------
			//	Particle relaxation starts here.
			//----------------------------------------------------------------------
			//random_wall_boundary_particles.parallel_exec(0.25);
			random_water_particles.parallel_exec(0.25);
			//random_vapor_particles.parallel_exec(0.25);
			//wall_boundary_relaxation_step_inner.surface_bounding_.parallel_exec();
			water_relaxation_step_inner.surface_bounding_.parallel_exec();
			//vapor_relaxation_step_inner.surface_bounding_.parallel_exec();
			update_smoothing_length_ratio_water.parallel_exec();
			water_block.updateCellLinkedListWithParticleSort(100);
			//write_wall_boundary_to_vtp.writeToFile(0);
			write_water_to_vtp.writeToFile(0);
			//write_vapor_to_vtp.writeToFile(0);
			//----------------------------------------------------------------------
			//	Relax particles of the insert body.
			//----------------------------------------------------------------------
			int ite_p = 0;
			while (ite_p < 1000)
			{
				//wall_boundary_relaxation_step_inner.parallel_exec();
				update_smoothing_length_ratio_water.parallel_exec();
				water_relaxation_step_inner.parallel_exec();
				//update_smoothing_length_ratio_vapor.parallel_exec();
				//vapor_relaxation_step_inner.parallel_exec();
				ite_p += 1;
				if (ite_p % 200 == 0)
				{
					std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the tank N = " << ite_p << "\n";
					write_water_to_vtp.writeToFile(ite_p);
					//write_vapor_to_vtp.writeToFile(ite_p);
					//write_wall_boundary_to_vtp.writeToFile(ite_p);
				}
			}
			std::cout << "The physics relaxation process of finish !" << std::endl;
			/** Output results. */
			//write_wall_boundary_particle_reload_files.writeToFile(0);
			write_water_particle_reload_files.writeToFile(0);
			//write_vapor_particle_reload_files.writeToFile(0);
			return 0;
		}

	    //ContactRelation water_wall_contact(water_block, { &wall_boundary });
		//AdaptiveContactRelation vapor_wall_contact(vapor_block, { &wall_boundary });
		ComplexRelation water_wall_complex(water_inner, /*{&vapor_block}*/{ &wall_boundary });
		//ComplexRelation vapor_water_complex(vapor_inner, /*{&water_block}*/water_wall_contact);


		//----------------------------------------------------------------------
		//	Define the numerical methods used in the simulation.
		//	Note that there may be data dependence on the sequence of constructions.
		//----------------------------------------------------------------------
		/** Initialize particle acceleration. */
		SimpleDynamics<NormalDirectionFromShapeAndOp> inner_normal_direction(wall_boundary, "InnerWall");
		SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, -gravity_g));
		SimpleDynamics<TimeStepInitialization> initialize_a_water_step(water_block,makeShared<VariableGravity>());
		//SimpleDynamics<TimeStepInitialization> initialize_a_vapor_step(vapor_block, makeShared<VariableGravity>());

		/** Evaluation of density by summation approach. */
		InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex>
			update_water_density_by_summation(water_wall_complex);
		//InteractionWithUpdate<fluid_dynamics::DensitySummationComplexAdaptive>
		//	update_vapor_density_by_summation(vapor_wall_contact,vapor_water_complex);
	//	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplexAdaptive>
		//	vapor_transport_correction(vapor_wall_contact,vapor_water_complex);
		/** Time step size without considering sound wave speed. */
		InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex>
			free_surface_indicator(water_wall_complex);
		InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplexAdaptive>
			water_transport_correction(water_wall_complex);
		ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_water_advection_time_step_size(water_block, U_max);
		//ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_vapor_advection_time_step_size(vapor_block, U_max);
		/** Time step size with considering sound wave speed. */
		ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_water_time_step_size(water_block);
		//	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_vapor_time_step_size(vapor_block);
			/** Pressure relaxation for water by using position verlet time stepping. */
		Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall>
			water_pressure_relaxation(water_wall_complex);
		Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall>
			water_density_relaxation(water_wall_complex);
		/** Extend Pressure relaxation is used for air. */
		//Dynamics1Level<fluid_dynamics::ExtendMultiPhaseIntegration1stHalfRiemannWithWall>
			//vapor_pressure_relaxation(vapor_wall_contact, vapor_water_complex,  2.0);
	//	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
		//	vapor_density_relaxation(vapor_wall_contact,vapor_water_complex);

			/** WaveProbes. */
		BodyRegionByCell wave_probe_buffer_no_4(water_block, makeShared<MultiPolygonShape>(createWaveProbeShape4(), "WaveProbe_04"));
		ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight, BodyRegionByCell>>
			wave_probe_4(io_environment, wave_probe_buffer_no_4);
		//----------------------------------------------------------------------
		//	Define the methods for I/O operations, observations
		//	and regression tests of the simulation.
		//----------------------------------------------------------------------
		BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
		RestartIO restart_io(io_environment, sph_system.real_bodies_);

		//----------------------------------------------------------------------
		//	Prepare the simulation with cell linked list, configuration
		//	and case specified initial condition if necessary.
		//----------------------------------------------------------------------
		sph_system.initializeSystemCellLinkedLists();
		sph_system.initializeSystemConfigurations();
		inner_normal_direction.parallel_exec();
		//----------------------------------------------------------------------
		//	Load restart file if necessary.
		//----------------------------------------------------------------------
		if (sph_system.restart_step_ != 0)
		{
			GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.restart_step_);
			water_block.updateCellLinkedList();
			//vapor_block.updateCellLinkedList();
			water_wall_complex.updateConfiguration();
			//vapor_water_complex.updateConfiguration();
		//	water_wall_contact.updateConfiguration();
			//vapor_wall_contact.updateConfiguration();

			//fluid_observer_contact.updateConfiguration();
		}
		//----------------------------------------------------------------------
		//	Setup for time-stepping control
		//----------------------------------------------------------------------
		size_t number_of_iterations = sph_system.restart_step_;
		int screen_output_interval = 100;
		int observation_sample_interval = screen_output_interval * 2;
		int restart_output_interval = screen_output_interval * 10;
		Real end_time = 3.0;
		Real output_interval = 0.05;
		Real dt = 0.0;
		//----------------------------------------------------------------------
		//	Statistics for CPU time
		//----------------------------------------------------------------------
		tick_count t1 = tick_count::now();
		tick_count::interval_t interval;
		tick_count::interval_t interval_computing_time_step;
		tick_count::interval_t interval_computing_fluid_pressure_relaxation;
		tick_count::interval_t interval_updating_configuration;
		tick_count time_instance;
		//----------------------------------------------------------------------
		//	First output before the main loop.
		//----------------------------------------------------------------------

		//thermo_water_initial_condition.parallel_exec();
		//thermo_vapor_initial_condition.parallel_exec();
		body_states_recording.writeToFile(0);
		wave_probe_4.writeToFile(0);
		//----------------------------------------------------------------------
		//	Main loop starts here.
		//----------------------------------------------------------------------

		while (GlobalStaticVariables::physical_time_ < end_time)
		{
			Real integration_time = 0.0;
			/** Integrate time (loop) until the next output time. */
			while (integration_time < output_interval)
			{
				/** Acceleration due to viscous force and gravity. */
				time_instance = tick_count::now();
				initialize_a_water_step.parallel_exec();
				//initialize_a_vapor_step.parallel_exec();
				free_surface_indicator.parallel_exec();
				Real Dt = get_water_advection_time_step_size.parallel_exec();
				//Real Dt_a = get_vapor_advection_time_step_size.parallel_exec();
				//Real Dt = SMIN(Dt_f, Dt_a);

				update_water_density_by_summation.parallel_exec();
				//	update_vapor_density_by_summation.parallel_exec();
				water_transport_correction.parallel_exec();
				//	vapor_transport_correction.parallel_exec(Dt);

				interval_computing_time_step += tick_count::now() - time_instance;

				/** Dynamics including pressure relaxation. */
				time_instance = tick_count::now();
				Real relaxation_time = 0.0;
				while (relaxation_time < Dt)
				{
					Real dt_f = get_water_time_step_size.parallel_exec();
					//	Real dt_a = get_vapor_time_step_size.parallel_exec();

					dt = SMIN(dt_f, Dt);

					water_pressure_relaxation.parallel_exec(dt);
					//	vapor_pressure_relaxation.parallel_exec(dt);

					water_density_relaxation.parallel_exec(dt);
					//	vapor_density_relaxation.parallel_exec(dt);

					relaxation_time += dt;
					integration_time += dt;
					GlobalStaticVariables::physical_time_ += dt;
				}
				interval_computing_fluid_pressure_relaxation += tick_count::now() - time_instance;

				if (number_of_iterations % screen_output_interval == 0)
				{
					std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						<< GlobalStaticVariables::physical_time_
						<< "	Dt = " << Dt << "	dt = " << dt << "\n";

					if (number_of_iterations % restart_output_interval == 0)
						restart_io.writeToFile(number_of_iterations);
				}
				number_of_iterations++;

				/** Update cell linked list and configuration. */
				time_instance = tick_count::now();

				water_block.updateCellLinkedListWithParticleSort(100);
				water_wall_complex.updateConfiguration();
				//water_wall_contact.updateConfiguration();

				//	vapor_block.updateCellLinkedListWithParticleSort(100);
				//	vapor_water_complex.updateConfiguration();
				//	vapor_wall_contact.updateConfiguration();

				interval_updating_configuration += tick_count::now() - time_instance;
			}

			tick_count t2 = tick_count::now();
			body_states_recording.writeToFile();
			wave_probe_4.writeToFile();
			tick_count t3 = tick_count::now();
			interval += t3 - t2;
		}

		tick_count t4 = tick_count::now();

		tick_count::interval_t tt;
		tt = t4 - t1 - interval;
		std::cout << "Total wall time for computation: " << tt.seconds()
			<< " seconds." << std::endl;
		std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
			<< interval_computing_time_step.seconds() << "\n";
		std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
			<< interval_computing_fluid_pressure_relaxation.seconds() << "\n";
		std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
			<< interval_updating_configuration.seconds() << "\n";



		return 0;
};
