/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * --------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
 * and HU1527/12-1.															*
 *                                                                           *
 * Portions copyright (c) 2017-2020 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * --------------------------------------------------------------------------*/
/**
 * @file 	solid_particles.h
 * @brief 	This is the derived class of base particles.
 * @author	Xiangyu Hu and Chi Zhang
 */

#ifndef SOLID_PARTICLES_H
#define SOLID_PARTICLES_H

#include "elastic_solid.h"
#include "base_particles.h"
#include "base_particles.hpp"

#include "particle_generator_lattice.h"
namespace SPH
{
	//----------------------------------------------------------------------
	//		preclaimed classes
	//----------------------------------------------------------------------
	class Solid;
	class ElasticSolid;

	/**
	 * @class SolidParticles
	 * @brief A group of particles with solid body particle data.
	 */
	class SolidParticles : public BaseParticles
	{
	public:
		SolidParticles(SPHBody &sph_body, Solid *solid);
		virtual ~SolidParticles(){};

		StdLargeVec<Vecd> pos0_; /**< initial position */
		StdLargeVec<Vecd> n_;	 /**< normal direction */
		StdLargeVec<Vecd> n0_;	 /**< initial normal direction */
		StdLargeVec<Matd> B_;	 /**< configuration correction for linear reproducing */
		//----------------------------------------------------------------------
		//		for fluid-structure interaction (FSI)
		//----------------------------------------------------------------------
		StdLargeVec<Vecd> vel_ave_;			 /**<  fluid time-step averaged particle velocity */
		StdLargeVec<Vecd> acc_ave_;			 /**<  fluid time-step averaged particle acceleration */
		StdLargeVec<Vecd> force_from_fluid_; /**<  forces (including pressure and viscous) from fluid */

		/** Get the kernel gradient in weak form. */
		virtual Vecd getKernelGradient(size_t index_i, size_t index_j, Real dW_ij, Vecd &e_ij) override;

		virtual void initializeOtherVariables() override;
		virtual SolidParticles *ThisObjectPtr() override { return this; };
	};

	/**
	 * @class ElasticSolidParticles
	 * @brief A group of particles with elastic body particle data.
	 */
	class ElasticSolidParticles : public SolidParticles
	{
	public:
		ElasticSolidParticles(SPHBody &sph_body, ElasticSolid *elastic_solid);
		virtual ~ElasticSolidParticles(){};

		StdLargeVec<Matd> F_;		   /**<  deformation tensor */
		StdLargeVec<Matd> dF_dt_;	   /**<  deformation tensor change rate */

		// STRAIN
		Matd getGreenLagrangeStrain(size_t particle_i);
		/**< Computing principal strain - returns the principal strains in descending order (starting from the largest) */
		Vecd getPrincipalStrains(size_t particle_i);
		/**< Computing von Mises equivalent strain from a static (constant) formulation. */
		Real getVonMisesStrain(size_t particle_i);
		/**< Computing von Mises equivalent strain from a "dynamic" formulation. This depends on the Poisson's ratio (from commercial FEM software Help). */
		Real getVonMisesStrainDynamic(size_t particle_i, Real poisson);

		/**< Computing von Mises strain for all particles. - "static" or "dynamic"*/
		StdLargeVec<Real> getVonMisesStrainVector(std::string strain_measure = "static");
		/**< Computing maximum von Mises strain from all particles. - "static" or "dynamic" */
		Real getVonMisesStrainMax(std::string strain_measure = "static");
		Real getPrincipalStrainMax();

		// STRESS
		Matd getStressCauchy(size_t particle_i);
		Matd getStressPK2(size_t particle_i);
		/**< Computing principal_stresses - returns the principal stresses in descending order (starting from the largest) */
		Vecd getPrincipalStresses(size_t particle_i);
		/**< Computing von_Mises_stress - "Cauchy" or "PK2" decided based on the stress_measure_ */
		Real getVonMisesStress(size_t particle_i);

		/**< Computing von Mises stress for all particles. - "Cauchy" or "PK2" decided based on the stress_measure_ */
		StdLargeVec<Real> getVonMisesStressVector();
		/**< Computing maximum von Mises stress from all particles. - "Cauchy" or "PK2" decided based on the stress_measure_ */
		Real getVonMisesStressMax();
		Real getPrincipalStressMax();

		/**< Computing displacement. */
		Vecd displacement(size_t particle_i);
		StdLargeVec<Vecd> getDisplacement();
		Real getMaxDisplacement();

		/**< Computing normal vector. */
		Vecd normal(size_t particle_i);
		StdLargeVec<Vecd> getNormal();

		/** relevant stress measure */
		std::string stress_measure_;

		virtual void initializeOtherVariables() override;
		virtual ElasticSolidParticles *ThisObjectPtr() override { return this; };

	protected:
		ElasticSolid *elastic_solid_;
	};

	/**
	 * @class ShellParticles
	 * @brief A group of particles with shell particle data.
	 */
	class ShellParticles : public ElasticSolidParticles
	{
	public:
		ShellParticles(SPHBody &sph_body, ElasticSolid *elastic_solid);
		virtual ~ShellParticles(){};

		Real thickness_ref_;
		StdLargeVec<Matd> transformation_matrix_; /**< initial transformation matrix from global to local coordinates */
		StdLargeVec<Real> thickness_;			  /**< shell thickness */
		//----------------------------------------------------------------------
		//	extra generalized coordinates in global coordinate
		//----------------------------------------------------------------------
		StdLargeVec<Vecd> pseudo_n_;	  /**< current pseudo-normal vector */
		StdLargeVec<Vecd> dpseudo_n_dt_;  /**< pseudo-normal vector change rate */
		StdLargeVec<Vecd> dpseudo_n_d2t_; /**< pseudo-normal vector second order time derivation */
		//----------------------------------------------------------------------
		//	extra generalized coordinate and velocity in local coordinate
		//----------------------------------------------------------------------
		StdLargeVec<Vecd> rotation_;		/**< rotation angle of the initial normal respective to each axis */
		StdLargeVec<Vecd> angular_vel_;		/**< angular velocity respective to each axis */
		StdLargeVec<Vecd> dangular_vel_dt_; /**< angular acceleration of respective to each axis*/
		//----------------------------------------------------------------------
		//	extra deformation and deformation rate in local coordinate
		//----------------------------------------------------------------------
		StdLargeVec<Matd> F_bending_;	  /**< bending deformation gradient	*/
		StdLargeVec<Matd> dF_bending_dt_; /**< bending deformation gradient change rate	*/
		//----------------------------------------------------------------------
		//	extra stress for pair interaction in global coordinate
		//----------------------------------------------------------------------
		StdLargeVec<Vecd> global_shear_stress_; /**< global shear stress */
		StdLargeVec<Matd> global_stress_;		/**<  global stress for pair interaction */
		StdLargeVec<Matd> global_moment_;		/**<  global bending moment for pair interaction */

		virtual void initializeOtherVariables() override;
		virtual ShellParticles *ThisObjectPtr() override { return this; };
	};
}
#endif // SOLID_PARTICLES_H
