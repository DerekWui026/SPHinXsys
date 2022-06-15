/**
 * @file 	elastic_dynamics.cpp
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "elastic_dynamics.h"
#include "general_dynamics.h"

#include <numeric>

using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//=================================================================================================//
		AcousticTimeStepSize::AcousticTimeStepSize(SPHBody &sph_body, Real CFL)
			: LocalParticleDynamicsReduce<Real, ReduceMin>(sph_body, DBL_MAX),
			  ElasticSolidDataSimple(sph_body), CFL_(CFL),
			  vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			  smoothing_length_(body_->sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//
		Real AcousticTimeStepSize::reduceRange(const IndexRange &particle_range, Real dt)
		{
			Real temp = reference_;
			for (size_t index_i = particle_range.begin(); index_i != particle_range.end(); ++index_i)
			{
				Real sound_speed = material_->ReferenceSoundSpeed();
				temp = operation_(temp, SMIN(sqrt(smoothing_length_ / (dvel_dt_[index_i].norm() + TinyReal)),
											 smoothing_length_ / (sound_speed + vel_n_[index_i].norm())));
			}
			return temp;
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::outputResult(Real reduced_value)
		{
			// since the particle does not change its configuration in pressure relaxation step
			// I chose a time-step size according to Eulerian method
			return CFL_ * reduced_value;
		}
		//=================================================================================================//
		ElasticDynamicsInitialCondition::
			ElasticDynamicsInitialCondition(SolidBody &solid_body)
			: LocalParticleDynamics(solid_body),
			  ElasticSolidDataSimple(solid_body),
			  pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_)
		{
		}
		//=================================================================================================//
		UpdateElasticNormalDirection::
			UpdateElasticNormalDirection(SolidBody &solid_body)
			: OldParticleDynamicsSimple(solid_body),
			  ElasticSolidDataSimple(solid_body),
			  n_(particles_->n_), n_0_(particles_->n_0_), F_(particles_->F_)
		{
		}
		//=================================================================================================//
		DeformationGradientTensorBySummation::
			DeformationGradientTensorBySummation(BaseBodyRelationInner &inner_relation)
			: OldInteractionDynamics(*inner_relation.sph_body_),
			  ElasticSolidDataInner(inner_relation),
			  Vol_(particles_->Vol_), pos_n_(particles_->pos_n_),
			  B_(particles_->B_), F_(particles_->F_)
		{
		}
		//=================================================================================================//
		void DeformationGradientTensorBySummation::Interaction(size_t index_i, Real dt)
		{
			Vecd &pos_n_i = pos_n_[index_i];

			Matd deformation(0.0);
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				deformation -= Vol_[index_j] * SimTK::outer((pos_n_i - pos_n_[index_j]), gradW_ij);
			}

			F_[index_i] = B_[index_i] * deformation;
		}
		//=================================================================================================//
		BaseElasticRelaxation::
			BaseElasticRelaxation(BaseBodyRelationInner &inner_relation)
			: LocalParticleDynamics(*inner_relation.sph_body_),
			  ElasticSolidDataInner(inner_relation), Vol_(particles_->Vol_),
			  rho_n_(particles_->rho_n_), mass_(particles_->mass_),
			  pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			  B_(particles_->B_), F_(particles_->F_), dF_dt_(particles_->dF_dt_) {}
		//=================================================================================================//
		StressRelaxationFirstHalf::
			StressRelaxationFirstHalf(BaseBodyRelationInner &inner_relation)
			: BaseElasticRelaxation(inner_relation),
			  dvel_dt_prior_(particles_->dvel_dt_prior_), force_from_fluid_(particles_->force_from_fluid_),
			  stress_PK1_(particles_->stress_PK1_), rho0_(material_->ReferenceDensity()), inv_rho0_(1.0 / rho0_),
			  smoothing_length_(body_->sph_adaptation_->ReferenceSmoothingLength()),
			  numerical_dissipation_factor_(0.25) {}
		//=================================================================================================//
		void StressRelaxationFirstHalf::initializeRange(const IndexRange &particle_range, Real dt)
		{
			for (size_t index_i = particle_range.begin(); index_i != particle_range.end(); ++index_i)
			{
				pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
			}

			for (size_t index_i = particle_range.begin(); index_i != particle_range.end(); ++index_i)
			{
				F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			}

			for (size_t index_i = particle_range.begin(); index_i != particle_range.end(); ++index_i)
			{
				rho_n_[index_i] = rho0_ / det(F_[index_i]);
			}

			for (size_t index_i = particle_range.begin(); index_i != particle_range.end(); ++index_i)
			{
				// obtain the first Piola-Kirchhoff stress from the second Piola-Kirchhoff stress
				// it seems using reproducing correction here increases convergence rate near the free surface
				stress_PK1_[index_i] = F_[index_i] * material_->ConstitutiveRelation(F_[index_i], index_i) * B_[index_i];
			}
		}
		//=================================================================================================//
		void StressRelaxationFirstHalf::interaction(size_t index_i, Real dt)
		{
			// including gravity and force from fluid
			Vecd acceleration = dvel_dt_prior_[index_i] + force_from_fluid_[index_i] / mass_[index_i];
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd e_ij = inner_neighborhood.e_ij_[n];
				Real r_ij = inner_neighborhood.r_ij_[n];
				Real dim_r_ij_1 = Dimensions / r_ij;
				Vecd pos_jump = pos_n_[index_i] - pos_n_[index_j];
				Vecd vel_jump = vel_n_[index_i] - vel_n_[index_j];
				Real strain_rate = SimTK::dot(pos_jump, vel_jump) * dim_r_ij_1 * dim_r_ij_1;
				Real weight = inner_neighborhood.W_ij_[n] * inv_W0_;
				Matd numerical_stress_ij =
					0.5 * (F_[index_i] + F_[index_j]) * material_->PairNumericalDamping(strain_rate, smoothing_length_);
				acceleration += (stress_PK1_[index_i] + stress_PK1_[index_j] +
								 numerical_dissipation_factor_ * weight * numerical_stress_ij) *
								inner_neighborhood.dW_ij_[n] * e_ij * Vol_[index_j] * inv_rho0_;
			}

			dvel_dt_[index_i] = acceleration;
		}
		//=================================================================================================//
		void StressRelaxationFirstHalf::updateRange(const IndexRange &particle_range, Real dt)
		{
			for (size_t index_i = particle_range.begin(); index_i != particle_range.end(); ++index_i)
			{
				vel_n_[index_i] += dvel_dt_[index_i] * dt;
			}
		}
		//=================================================================================================//
		KirchhoffParticleStressRelaxationFirstHalf::
			KirchhoffParticleStressRelaxationFirstHalf(BaseBodyRelationInner &inner_relation)
			: StressRelaxationFirstHalf(inner_relation){};
		//=================================================================================================//
		void KirchhoffParticleStressRelaxationFirstHalf::
			initializeRange(const IndexRange &particle_range, Real dt)
		{
			for (size_t index_i = particle_range.begin(); index_i != particle_range.end(); ++index_i)
			{
				pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
			}

			for (size_t index_i = particle_range.begin(); index_i != particle_range.end(); ++index_i)
			{
				F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			}

			for (size_t index_i = particle_range.begin(); index_i != particle_range.end(); ++index_i)
			{
				Real J = det(F_[index_i]);
				Real one_over_J = 1.0 / J;
				rho_n_[index_i] = rho0_ * one_over_J;
				Real J_to_minus_2_over_dimension = pow(one_over_J, 2.0 * one_over_dimensions_);
				Matd normalized_b = (F_[index_i] * ~F_[index_i]) * J_to_minus_2_over_dimension;
				Matd deviatoric_b = normalized_b - Matd(1.0) * normalized_b.trace() * one_over_dimensions_;
				Matd inverse_F_T = ~SimTK::inverse(F_[index_i]);
				// obtain the first Piola-Kirchhoff stress from the Kirchhoff stress
				// it seems using reproducing correction here increases convergence rate
				// near the free surface however, this correction is not used for the numerical dissipation
				stress_PK1_[index_i] = (Matd(1.0) * material_->VolumetricKirchhoff(J) +
										material_->DeviatoricKirchhoff(deviatoric_b)) *
									   inverse_F_T * B_[index_i];
			}
		}
		//=================================================================================================//
		KirchhoffStressRelaxationFirstHalf::
			KirchhoffStressRelaxationFirstHalf(BaseBodyRelationInner &inner_relation)
			: StressRelaxationFirstHalf(inner_relation)
		{
			particles_->registerAVariable(J_to_minus_2_over_dimension_, "DeterminantTerm");
			particles_->registerAVariable(stress_on_particle_, "StressOnParticle");
		};
		//=================================================================================================//
		void KirchhoffStressRelaxationFirstHalf::initializeRange(const IndexRange &particle_range, Real dt)
		{
			for (size_t index_i = particle_range.begin(); index_i != particle_range.end(); ++index_i)
			{
				pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
			}

			for (size_t index_i = particle_range.begin(); index_i != particle_range.end(); ++index_i)
			{
				F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			}

			for (size_t index_i = particle_range.begin(); index_i != particle_range.end(); ++index_i)
			{
				Real J = det(F_[index_i]);
				Real one_over_J = 1.0 / J;
				rho_n_[index_i] = rho0_ * one_over_J;
				J_to_minus_2_over_dimension_[index_i] = pow(one_over_J * one_over_J, one_over_dimensions_);
				Matd inverse_F_T = ~SimTK::inverse(F_[index_i]);
				stress_on_particle_[index_i] =
					inverse_F_T * (material_->VolumetricKirchhoff(J) -
								   correction_factor_ * material_->ShearModulus() * J_to_minus_2_over_dimension_[index_i] *
									   (F_[index_i] * ~F_[index_i]).trace() * one_over_dimensions_) +
					material_->NumericalDampingLeftCauchy(F_[index_i], dF_dt_[index_i], smoothing_length_, index_i) * inverse_F_T;
			}

			for (size_t index_i = particle_range.begin(); index_i != particle_range.end(); ++index_i)
			{
				stress_PK1_[index_i] = F_[index_i] * material_->ConstitutiveRelation(F_[index_i], index_i);
			}
		}
		//=================================================================================================//
		void KirchhoffStressRelaxationFirstHalf::interaction(size_t index_i, Real dt)
		{
			// including gravity and force from fluid
			Vecd acceleration = dvel_dt_prior_[index_i] + force_from_fluid_[index_i] / mass_[index_i];
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd shear_force_ij = correction_factor_ * material_->ShearModulus() *
									  (J_to_minus_2_over_dimension_[index_i] + J_to_minus_2_over_dimension_[index_j]) *
									  (pos_n_[index_i] - pos_n_[index_j]) / inner_neighborhood.r_ij_[n];
				acceleration += ((stress_on_particle_[index_i] + stress_on_particle_[index_j]) * inner_neighborhood.e_ij_[n] + shear_force_ij) *
								inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inv_rho0_;
			}
			dvel_dt_[index_i] = acceleration;
		}
		//=================================================================================================//
		void StressRelaxationSecondHalf::initializeRange(const IndexRange &particle_range, Real dt)
		{
			for (size_t index_i = particle_range.begin(); index_i != particle_range.end(); ++index_i)
			{
				pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
			}
		}
		//=================================================================================================//
		void StressRelaxationSecondHalf::interaction(size_t index_i, Real dt)
		{
			const Vecd &vel_n_i = vel_n_[index_i];

			dF_dt_[index_i] = Matd(0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				dF_dt_[index_i] -= Vol_[index_j] * SimTK::outer((vel_n_i - vel_n_[index_j]), gradW_ij);
			}
		}
		//=================================================================================================//
		void StressRelaxationSecondHalf::updateRange(const IndexRange &particle_range, Real dt)
		{
			for (size_t index_i = particle_range.begin(); index_i != particle_range.end(); ++index_i)
			{
				dF_dt_[index_i] *= B_[index_i];
			}

			for (size_t index_i = particle_range.begin(); index_i != particle_range.end(); ++index_i)
			{
				F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			}
		}
		//=================================================================================================//
	}
}
