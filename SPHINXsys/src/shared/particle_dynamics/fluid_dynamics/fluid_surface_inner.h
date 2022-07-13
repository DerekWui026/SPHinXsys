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
* @file 	fluid_surface_inner.h
* @brief 	Here, we define the algorithm classes for fluid surfaces. 
* @details 	Fluid indicators are mainly used here to classify different region in a fluid.   
* @author	Chi Zhang and Xiangyu Hu
*/

#ifndef FLUID_SURFACE_INNER_H
#define FLUID_SURFACE_INNER_H

#include "fluid_dynamics_inner.hpp"

namespace SPH
{
	namespace fluid_dynamics
	{
		/**
		* @class FreeSurfaceIndicationInner
		* @brief  indicate the particles near the free surface of a fluid body.
		* Note that, SPHinXsys does not require this function for simulating general free surface flow problems.
		* However, some other applications may use this function, such as transport velocity formulation, 
		* for masking some function which is only applicable for the bulk of the fluid body.
		*/
		class FreeSurfaceIndicationInner
			: public InteractionDynamicsWithUpdate,
			  public FluidDataInner
		{
		public:
			explicit FreeSurfaceIndicationInner(BaseBodyRelationInner &inner_relation, Real threshold = 0.75);
			virtual ~FreeSurfaceIndicationInner(){};

		protected:
			Real threshold_by_dimensions_;
			StdLargeVec<Real> &Vol_;
			StdLargeVec<int> &surface_indicator_;
			StdLargeVec<Real> pos_div_;
			Real smoothing_length_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
        * @class SpatialTemporalFreeSurfaceIdentification
        * @brief using the spatial-temporal method to indicate the surface particles to avoid mis-judgement.
        */
		template <class FreeSurfaceIdentification>
		class SpatialTemporalFreeSurfaceIdentification : public FreeSurfaceIdentification
		{
		public:
			template <typename... ConstructorArgs>
			explicit SpatialTemporalFreeSurfaceIdentification(ConstructorArgs &&...args);
			virtual ~SpatialTemporalFreeSurfaceIdentification(){};

		protected:
			StdLargeVec<int> previous_surface_indicator_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;

			void checkNearPreviousFreeSurface(size_t index_i);
		};
		using SpatialTemporalFreeSurfaceIdentificationInner =
			SpatialTemporalFreeSurfaceIdentification<FreeSurfaceIndicationInner>;

		/**
		 * @class DensitySummationFreeSurfaceInner
		 * @brief computing density by summation with a re-normalization for free surface flows
		 */
		class DensitySummationFreeSurfaceInner : public DensitySummationInner
		{
		public:
			explicit DensitySummationFreeSurfaceInner(BaseBodyRelationInner &inner_relation)
				: DensitySummationInner(inner_relation){};
			virtual ~DensitySummationFreeSurfaceInner(){};

		protected:
			virtual Real ReinitializedDensity(Real rho_sum, Real rho_0, Real rho_n) override
			{
				return rho_sum + SMAX(0.0, (rho_n - rho_sum)) * rho_0 / rho_n;
			};
		};

		/**
		 * @class DensitySummationFreeStreamInner
		 * @brief The density is smoothed if the particle is near fluid surface.
		 */
		class DensitySummationFreeStreamInner : public DensitySummationFreeSurfaceInner
		{
		public:
			explicit DensitySummationFreeStreamInner(BaseBodyRelationInner &inner_relation)
				: DensitySummationFreeSurfaceInner(inner_relation),
				  surface_indicator_(*particles_->getVariableByName<int>("SurfaceIndicator")){};
			virtual ~DensitySummationFreeStreamInner(){};

		protected:
			StdLargeVec<int> &surface_indicator_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
			bool isNearSurface(size_t index_i);
		};

		/**
        * @class FreeStreamBoundaryVelocityCorrection
        * @brief this function is applied to freestream flows
		* @brief modify the velocity of free surface particles with far-field velocity
        */
		class FreeStreamBoundaryVelocityCorrection : public ParticleDynamicsSimple, public FluidDataInner
		{
		public:
			explicit FreeStreamBoundaryVelocityCorrection(BaseBodyRelationInner &inner_relation)
				: ParticleDynamicsSimple(*inner_relation.sph_body_),
				  FluidDataInner(inner_relation), u_ref_(1.0), t_ref_(2.0),
				  rho_ref_(material_->ReferenceDensity()), rho_sum(particles_->rho_sum_),
				  vel_(particles_->vel_), acc_(particles_->acc_),
				  surface_indicator_(*particles_->getVariableByName<int>("SurfaceIndicator")){};
			virtual ~FreeStreamBoundaryVelocityCorrection(){};

		protected:
			Real u_ref_, t_ref_, rho_ref_;
			StdLargeVec<Real> &rho_sum;
			StdLargeVec<Vecd> &vel_, &acc_;
			StdLargeVec<int> &surface_indicator_;

			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class FreeSurfaceProbeOnFluidBody
		 * @brief Probe the free surface profile for a fluid body part by reduced operation.
		 */
		class FreeSurfaceProbeOnFluidBody : public PartDynamicsByCellReduce<Real, ReduceMax>,
											public FluidDataSimple
		{
		public:
			FreeSurfaceProbeOnFluidBody(FluidBody &fluid_body, BodyPartByCell &body_part)
				: PartDynamicsByCellReduce<Real, ReduceMax>(fluid_body, body_part), FluidDataSimple(fluid_body),
				  pos_(particles_->pos_)
			{
				quantity_name_ = body_part.BodyPartName() + "FreeSurface";
				initial_reference_ = 0.0;
			}
			virtual ~FreeSurfaceProbeOnFluidBody(){};

		protected:
			StdLargeVec<Vecd> &pos_;
			virtual void SetupReduce() override{};
			virtual Real ReduceFunction(size_t index_i, Real dt = 0.0) override { return pos_[index_i][1]; };
		};
		/**
		 * @class ColorFunctionGradientInner
		 * @brief  indicate the particles near the interface of a fluid-fluid interaction and computing norm
		 */
		class ColorFunctionGradientInner : public InteractionDynamics, public FluidDataInner
		{
		public:
			explicit ColorFunctionGradientInner(BaseBodyRelationInner &inner_relation);
			virtual ~ColorFunctionGradientInner(){};

		protected:
			Real threshold_by_dimensions_;
			StdLargeVec<Real> &Vol_;
			StdLargeVec<int> &surface_indicator_;
			StdLargeVec<Vecd> color_grad_;
			StdLargeVec<Vecd> surface_norm_;
			StdLargeVec<Real> &pos_div_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ColorFunctionGradientInterpolationInner
		 * @brief  the viscous force induced acceleration
		 */
		class ColorFunctionGradientInterpolationInner : public InteractionDynamics, public FluidDataInner
		{
		public:
			explicit ColorFunctionGradientInterpolationInner(BaseBodyRelationInner &inner_relation);
			virtual ~ColorFunctionGradientInterpolationInner(){};

		protected:
			Real threshold_by_dimensions_;
			StdLargeVec<Real> &Vol_;
			StdLargeVec<int> &surface_indicator_;
			StdLargeVec<Vecd> &color_grad_;
			StdLargeVec<Vecd> &surface_norm_;
			StdLargeVec<Real> &pos_div_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class SurfaceTensionAccelerationInner
		 * @brief  the viscous force induced acceleration
		 */
		class SurfaceTensionAccelerationInner : public InteractionDynamics, public FluidDataInner
		{
		public:
			SurfaceTensionAccelerationInner(BaseBodyRelationInner &inner_relation, Real gamma);
			explicit SurfaceTensionAccelerationInner(BaseBodyRelationInner &inner_relation);
			virtual ~SurfaceTensionAccelerationInner(){};

		protected:
			Real gamma_;
			StdLargeVec<Real> &Vol_, &mass_;
			StdLargeVec<Vecd> &acc_prior_;
			StdLargeVec<int> &surface_indicator_;
			StdLargeVec<Vecd> &color_grad_;
			StdLargeVec<Vecd> &surface_norm_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
	}
}
#endif //FLUID_SURFACE_INNER_H