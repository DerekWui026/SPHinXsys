/**
 * @file 	complex_solid.hpp
 * @brief 	These are implementation of the classes for complex solids.
 * @author	Xiangyu Hu and Chi Zhang
 */
#pragma once

#include "complex_solid.h"

#include "base_particles.hpp"

using namespace std;

namespace SPH
{
	//=============================================================================================//
	template <class MuscleType>
	template <typename... ConstructorArgs>
	ActiveMuscle<MuscleType>::ActiveMuscle(ConstructorArgs &&...args)
		: MuscleType(std::forward<ConstructorArgs>(args)...)
	{
		MuscleType::material_type_name_ = "ActiveMuscle";
	}
	//=============================================================================================//
	template <class MuscleType>
	void ActiveMuscle<MuscleType>::initializeContractionStress()
	{
		this->base_particles_->registerAVariable(active_contraction_stress_, "ActiveContractionStress");
	}
	//=============================================================================================//
	template <class MuscleType>
	void ActiveMuscle<MuscleType>::assignBaseParticles(BaseParticles *base_particles)
	{
		MuscleType::assignBaseParticles(base_particles);
		initializeContractionStress();
	}
	//=============================================================================================//
	template <class MuscleType>
	Matd ActiveMuscle<MuscleType>::ConstitutiveRelation(Matd &deformation, size_t index_i)
	{
		return MuscleType::ConstitutiveRelation(deformation, index_i) +
			   active_contraction_stress_[index_i] * MuscleType::MuscleFiberDirection(index_i);
	}
	//=============================================================================================//
}
