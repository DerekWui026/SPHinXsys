/**
 * @file 	inelastic_dynamics.cpp
 * @author	Xiaojing Tang, Chi Zhang and Xiangyu Hu
 */

#include "inelastic_dynamics.h"

using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//=================================================================================================//
		PlasticStressRelaxationFirstHalf::
			PlasticStressRelaxationFirstHalf(BaseBodyRelationInner &inner_relation) :
			StressRelaxationFirstHalf(inner_relation),
			plastic_solid_(DynamicCast<PlasticSolid>(this, material_))
		{
			numerical_dissipation_factor_ = 0.5;
		}
		//=================================================================================================//
		void PlasticStressRelaxationFirstHalf::Initialization(size_t index_i, Real dt)
		{
			pos_[index_i] += vel_[index_i] * dt * 0.5;
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			rho_[index_i] = rho0_ / SimTK::det(F_[index_i]);

			// TODO: this formulation is temporary, to be revises. 
			stress_PK1_B_[index_i] = plastic_solid_->PlasticConstitutiveRelation(F_[index_i], index_i, dt);
		}
		//=================================================================================================//
	}
}
