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
* @file 	general_solid_dynamics.h
* @brief 	Here, we define the algorithm classes for solid dynamics. 
* @details 	We consider here a weakly compressible solids.   
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
*/

#ifndef GENERAL_SOLID_DYNAMICS_H
#define GENERAL_SOLID_DYNAMICS_H

#include "all_particle_dynamics.h"
#include "general_dynamics.h"
#include "base_kernel.h"
#include "all_body_relations.h"
#include "solid_body.h"
#include "solid_particles.h"
#include "elastic_solid.h"



namespace SPH
{
	namespace solid_dynamics
	{
		//----------------------------------------------------------------------
		//		for general solid dynamics
		//----------------------------------------------------------------------
		typedef DataDelegateSimple<SolidBody, SolidParticles, Solid> SolidDataSimple;
		typedef DataDelegateInner<SolidBody, SolidParticles, Solid> SolidDataInner;

		/**
		 * @class SolidDynamicsInitialCondition
		 * @brief  set initial condition for solid fluid body
		 * This is a abstract class to be override for case specific initial conditions.
		 */
		class SolidDynamicsInitialCondition : public ParticleDynamicsSimple, public SolidDataSimple
		{
		public:
			explicit SolidDynamicsInitialCondition(SolidBody &solid_body)
				: ParticleDynamicsSimple(solid_body), SolidDataSimple(solid_body){};
			virtual ~SolidDynamicsInitialCondition(){};
		};

		/**
		* @class CorrectConfiguration
		* @brief obtain the corrected initial configuration in strong form
		*/
		class CorrectConfiguration : public InteractionDynamics, public SolidDataInner
		{
		public:
			explicit CorrectConfiguration(BaseBodyRelationInner &inner_relation);
			virtual ~CorrectConfiguration(){};

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Matd> &B_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

	}
}
#endif //GENERAL_SOLID_DYNAMICS_H