/**
 * @file 	3d_elasticSolid_shell_collision.cpp
 * @brief 	This is a benchmark test of the 3D elastic solid->shell contact/impact formulations.
 * @details  We consider the collision of an elastic ball bouncing in a spherical shell box.
 * @author 	Massoud Rezavand, Virtonomy GmbH
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real resolution_ref = 0.05;							/**< reference resolution. */
Real thickness = resolution_ref * 1.;				/**< shell thickness. */
Real radius = 2.0;									/**< cylinder radius. */
Real half_height = 1.0;								/** Height of the cylinder. */
Real radius_mid_surface = radius + thickness / 2.0; /** Radius of the mid surface. */
int particle_number_mid_surface = int(2.0 * radius_mid_surface * Pi * 215.0 / 360.0 / resolution_ref);
int particle_number_height = 2 * int(half_height / resolution_ref);
int BWD = 1;						  /** Width of the boundary layer measured by number of particles. */
BoundingBox system_domain_bounds(Vec3d(-radius - thickness, -half_height - thickness, -radius - thickness),
								 Vec3d(radius + thickness, half_height + thickness, radius + thickness));
Real ball_radius = 0.5;
Real gravity_g = 1.0;
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;
Real Youngs_modulus = 2.0e4;
Real poisson = 0.45;
Real physical_viscosity = 1.0e6;
//----------------------------------------------------------------------
/** Define application dependent particle generator for thin structure. */
class CylinderParticleGenerator : public SurfaceParticleGenerator
{
public:
	explicit CylinderParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
	virtual void initializeGeometricVariables() override
	{
		// the cylinder and boundary
		for (int i = 0; i < particle_number_mid_surface + 2 * BWD; i++)
		{
			for (int j = 0; j < particle_number_height; j++)
			{
				Real x = radius_mid_surface * cos(162.5 / 180.0 * Pi + (i - BWD + 0.5) * 215.0 / 360.0 * 2 * Pi / (Real)particle_number_mid_surface);
				Real y = (j -  particle_number_height / 2) * resolution_ref + resolution_ref * 0.5;
				Real z = radius_mid_surface * sin(162.5 / 180.0 * Pi + (i - BWD + 0.5) * 215.0 / 360.0 * 2 * Pi / (Real)particle_number_mid_surface);
				initializePositionAndVolumetricMeasure(Vecd(x, y, z), resolution_ref * resolution_ref);
				Vec3d n_0 = Vec3d(x / radius_mid_surface, 0.0, z / radius_mid_surface);
				initializeSurfaceProperties(n_0, thickness);
			}
		}
	}
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	/** Tag for running particle relaxation for the initially body-fitted distribution */
	sph_system.run_particle_relaxation_ = false;
	/** Tag for starting with relaxed body-fitted particles distribution */
	sph_system.reload_particles_ = true;
	/** Tag for computation from restart files. 0: start with initial condition */
	sph_system.restart_step_ = 0;
	/** Handle command line arguments. */
	sph_system.handleCommandlineOptions(ac, av);
	/** I/O environment. */
	InOutput in_output(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	/** create a shell body. */
	SolidBody shell(sph_system, makeShared<DefaultShape>("shell"));
	shell.defineParticlesAndMaterial<ShellParticles, LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
	shell.generateParticles<CylinderParticleGenerator>();

	SolidBody ball(sph_system, makeShared<GeometricShapeBall>(Vec3d(radius / 2.0, 0.0, 0.0), ball_radius, "BallBody"));
	ball.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
	if (!sph_system.run_particle_relaxation_ && sph_system.reload_particles_)
	{
		ball.generateParticles<ParticleGeneratorReload>(in_output, ball.getBodyName());
	}
	else
	{
		ball.defineBodyLevelSetShape()->writeLevelSet(ball);
		ball.generateParticles<ParticleGeneratorLattice>();
	}
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.run_particle_relaxation_)
	{
		//----------------------------------------------------------------------
		//	Define body relation map used for particle relaxation.
		//----------------------------------------------------------------------
		BodyRelationInner ball_inner(ball);
		//----------------------------------------------------------------------
		//	Define the methods for particle relaxation for ball.
		//----------------------------------------------------------------------
		RandomizeParticlePosition ball_random_particles(ball);
		relax_dynamics::RelaxationStepInner ball_relaxation_step_inner(ball_inner);
		//----------------------------------------------------------------------
		//	Output for particle relaxation.
		//----------------------------------------------------------------------
		BodyStatesRecordingToVtp write_relaxed_particles(in_output, sph_system.real_bodies_);
		ReloadParticleIO write_particle_reload_files(in_output, {&ball});
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		ball_random_particles.parallel_exec(0.25);
		write_relaxed_particles.writeToFile(0);
		//----------------------------------------------------------------------
		//	From here iteration for particle relaxation begins.
		//----------------------------------------------------------------------
		int ite = 0;
		int relax_step = 1000;
		while (ite < relax_step)
		{
			ball_relaxation_step_inner.parallel_exec();
			ite += 1;
			if (ite % 100 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
				write_relaxed_particles.writeToFile(ite);
			}
		}
		std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
		write_particle_reload_files.writeToFile(0);
		return 0;
	}
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner ball_inner(ball);
	SolidBodyRelationContact ball_contact(ball, {&shell});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	Gravity gravity(Vec3d(0.0, 0.0, -gravity_g));
	TimeStepInitialization ball_initialize_timestep(ball, gravity);
	solid_dynamics::CorrectConfiguration ball_corrected_configuration(ball_inner);
	solid_dynamics::AcousticTimeStepSize ball_get_time_step_size(ball, 0.45);
	/** stress relaxation for the balls. */
	solid_dynamics::KirchhoffStressRelaxationFirstHalf ball_stress_relaxation_first_half(ball_inner);
	solid_dynamics::StressRelaxationSecondHalf ball_stress_relaxation_second_half(ball_inner);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::ShellContactDensity ball_update_contact_density(ball_contact);
	solid_dynamics::ContactForceFromWall ball_compute_solid_contact_forces(ball_contact);
	DampingWithRandomChoice<solid_dynamics::PairwiseFrictionFromWall> ball_friction(0.1, ball_contact, physical_viscosity);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp body_states_recording(in_output, sph_system.real_bodies_);
	BodyStatesRecordingToVtp write_ball_state(in_output, {ball});
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	ball_corrected_configuration.parallel_exec();
	/** Initial states output. */
	body_states_recording.writeToFile(0);
	/** Main loop. */
	int ite = 0;
	Real T0 = 10.0;
	Real End_Time = T0;
	Real D_Time = 0.01 * T0;
	Real Dt = 0.1 * D_Time;
	Real dt = 0.0;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		while (integration_time < D_Time)
		{
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				ball_initialize_timestep.parallel_exec();
				if (ite % 100 == 0)
				{
					std::cout << "N=" << ite << " Time: "
							  << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
				}
				ball_update_contact_density.parallel_exec();
				ball_compute_solid_contact_forces.parallel_exec();
				ball_stress_relaxation_first_half.parallel_exec(dt);
				ball_friction.parallel_exec(dt);
				ball_stress_relaxation_second_half.parallel_exec(dt);

				ball.updateCellLinkedList();
				ball_contact.updateConfiguration();

				ite++;
				Real dt_free = ball_get_time_step_size.parallel_exec();
				dt = dt_free;
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
		}
		tick_count t2 = tick_count::now();
		write_ball_state.writeToFile(ite);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	return 0;
}
