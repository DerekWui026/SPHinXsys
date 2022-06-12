/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1                *
 * and HU1527/12-1.                                                             *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
* @file 	particle_dynamics_algorithms.h
* @brief 	This is the classes for algorithms particle dynamics.
* @detail	Generally, there are four types dynamics. One is without particle interaction.
*			One is with particle interaction within a body. One is with particle interaction
*			between a center body and other contacted bodies.
* 			Still another is the combination of the last two.
*			For the first dynamics, there is also reduce dynamics
*			which carries reduced operations through the particles of the body.

* @author	Chi ZHang and Xiangyu Hu
*/

#ifndef PARTICLE_DYNAMICS_ALGORITHMS_H
#define PARTICLE_DYNAMICS_ALGORITHMS_H

#include "base_particle_dynamics.hpp"

namespace SPH
{
	/**
	 * @class DynamicsParameter
	 * @brief Obtain a parameter without looping particles
	 */
	template <class LocalDynamicsType>
	class DynamicsParameter : public LocalDynamicsType,
							  public BaseParticleDynamics<
								  typename LocalDynamicsType::DynamicsParamterType>
	{
		using DynamicsParamterType = typename LocalDynamicsType::DynamicsParamterType;

	public:
		template <typename... Args>
		explicit DynamicsParameter(Args &&...args)
			: LocalDynamicsType(std::forward<Args>(args)...),
			  BaseParticleDynamics<typename LocalDynamicsType::DynamicsParamterType>(){};
		virtual ~DynamicsParameter(){};

		virtual DynamicsParamterType exec(Real dt = 0.0) override
		{
			LocalDynamicsType::setBodyUpdated();
			return LocalDynamicsType::setupDynamics(dt);
		};

		virtual DynamicsParamterType parallel_exec(Real dt = 0.0) override { return exec(dt); };
	};

	/**
	 * @class SimpleDynamics
	 * @brief Simple particle dynamics without considering particle interaction
	 */
	template <class LocalDynamicsType>
	class SimpleDynamics : public LocalDynamicsType,
						   public BaseSimpleDynamics<size_t, RangeFunctor>
	{
	public:
		template <typename... Args>
		explicit SimpleDynamics(Args &&...args)
			: LocalDynamicsType(std::forward<Args>(args)...),
			  BaseSimpleDynamics<size_t, RangeFunctor>(
				  LocalDynamicsType::body_->BodyRange(),
				  std::bind(&LocalDynamicsType::updateRange, this, _1, _2)){};
		virtual ~SimpleDynamics(){};

		virtual void runSetup(Real dt = 0.0) override
		{
			LocalDynamicsType::setBodyUpdated();
			LocalDynamicsType::setupDynamics(dt);
		};
	};

	/**
	 * @class InteractionDynamics
	 * @brief  This is the class for particle interaction with other particles
	 */
	template <class LocalDynamicsType>
	class InteractionDynamics : public LocalDynamicsType,
								public BaseInteractionDynamics<size_t>
	{
	public:
		template <typename... Args>
		explicit InteractionDynamics(Args &&...args)
			: LocalDynamicsType(std::forward<Args>(args)...),
			  BaseInteractionDynamics<size_t>(
				  LocalDynamicsType::body_->BodyRange(),
				  std::bind(&LocalDynamicsType::interaction, this, _1, _2)){};
		virtual ~InteractionDynamics(){};

		virtual void runSetup(Real dt = 0.0) override
		{
			LocalDynamicsType::setBodyUpdated();
			LocalDynamicsType::setupDynamics(dt);
		};
	};

	/**
	 * @class InteractionDynamics
	 * @brief This class includes an interaction and a update steps
	 */
	template <class LocalDynamicsType>
	class InteractionDynamicsAndUpdate : public LocalDynamicsType,
										 public BaseInteractionDynamicsAndUpdate<size_t, RangeFunctor>
	{
	public:
		template <typename... Args>
		explicit InteractionDynamicsAndUpdate(Args &&...args)
			: LocalDynamicsType(std::forward<Args>(args)...),
			  BaseInteractionDynamicsAndUpdate<size_t, RangeFunctor>(
				  LocalDynamicsType::body_->BodyRange(),
				  std::bind(&LocalDynamicsType::interaction, this, _1, _2),
				  std::bind(&LocalDynamicsType::updateRange, this, _1, _2)){};
		virtual ~InteractionDynamicsAndUpdate(){};

		virtual void runSetup(Real dt = 0.0) override
		{
			LocalDynamicsType::setBodyUpdated();
			LocalDynamicsType::setupDynamics(dt);
		};
	};

	/**
	 * @class InteractionDynamics1Level
	 * @brief This class includes an initialization, an interaction and a update steps
	 */
	template <class LocalDynamicsType>
	class InteractionDynamics1Level : public LocalDynamicsType,
									  public BaseInteractionDynamics1Level<size_t, RangeFunctor>
	{
	public:
		template <typename... Args>
		explicit InteractionDynamics1Level(Args &&...args)
			: LocalDynamicsType(std::forward<Args>(args)...),
			  BaseInteractionDynamics1Level<size_t, RangeFunctor>(
				  LocalDynamicsType::body_->BodyRange(),
				  std::bind(&LocalDynamicsType::initializeRange, this, _1, _2),
				  std::bind(&LocalDynamicsType::interaction, this, _1, _2),
				  std::bind(&LocalDynamicsType::updateRange, this, _1, _2)){};
		virtual ~InteractionDynamics1Level() override{};

		virtual void runSetup(Real dt = 0.0) override
		{
			LocalDynamicsType::setBodyUpdated();
			LocalDynamicsType::setupDynamics(dt);
		};
	};

	/**
	 * @class MultiStageDynamics
	 * @brief The algorithm inlcudes a initilization step,
	 * and serveral successive interaction steps.
	 */
	class MultiStageDynamics : public BaseParticleDynamics<void>
	{
		UniquePtrKeepers<BaseParticleDynamics<void> *> stages_keeper_;
		StdVec<BaseParticleDynamics<void> *> multi_stage_steps_;

	public:
		MultiStageDynamics() : BaseParticleDynamics<void>(){};
		virtual ~MultiStageDynamics(){};

		template <class ParticleDynamicsType, typename... Args>
		MultiStageDynamics *addParticleDyanmics(Args &&...args)
		{
			multi_stage_steps_.push_back(
				stages_keeper_.createPtr<ParticleDynamicsType>(std::forward<Args>(args)...));
			return this;
		};
		virtual void exec(Real dt = 0.0) override
		{
			for (size_t k = 0; k < multi_stage_steps_.size(); ++k)
				multi_stage_steps_[k]->exec(dt);
		};

		virtual void parallel_exec(Real dt = 0.0) override
		{
			for (size_t k = 0; k < multi_stage_steps_.size(); ++k)
				multi_stage_steps_[k]->parallel_exec(dt);
		};
	};

	/**
	 * @class MultipleInteractionDynamics1Level
	 * @brief Several body dynamic updating together so that they
	 * are not dependent on their sequence.
	 */
	class MultipleInteractionDynamics1Level : public BaseParticleDynamics<void>
	{
		StdVec<BaseInteractionDynamics1Level<size_t, RangeFunctor> *> bodies_dynamics_1level_;

	public:
		explicit MultipleInteractionDynamics1Level(
			StdVec<BaseInteractionDynamics1Level<size_t, RangeFunctor> *> bodies_dynamics_1level)
			: BaseParticleDynamics<void>(), bodies_dynamics_1level_(bodies_dynamics_1level){};
		virtual ~MultipleInteractionDynamics1Level(){};

		virtual void exec(Real dt = 0.0) override
		{
			for (size_t k = 0; k < bodies_dynamics_1level_.size(); ++k)
				bodies_dynamics_1level_[k]->runSetup(dt);
			for (size_t k = 0; k < bodies_dynamics_1level_.size(); ++k)
				bodies_dynamics_1level_[k]->runInitialization(dt);
			for (size_t k = 0; k < bodies_dynamics_1level_.size(); ++k)
				bodies_dynamics_1level_[k]->runInteraction(dt);
			for (size_t k = 0; k < bodies_dynamics_1level_.size(); ++k)
				bodies_dynamics_1level_[k]->runUpdate(dt);
		};

		virtual void parallel_exec(Real dt = 0.0) override
		{
			for (size_t k = 0; k < bodies_dynamics_1level_.size(); ++k)
				bodies_dynamics_1level_[k]->runSetup(dt);
			for (size_t k = 0; k < bodies_dynamics_1level_.size(); ++k)
				bodies_dynamics_1level_[k]->runInitialization_parallel(dt);
			for (size_t k = 0; k < bodies_dynamics_1level_.size(); ++k)
				bodies_dynamics_1level_[k]->runInteraction_parallel(dt);
			for (size_t k = 0; k < bodies_dynamics_1level_.size(); ++k)
				bodies_dynamics_1level_[k]->runUpdate_parallel(dt);
		};
	};

	/**
	 * @class SimpleDynamicsReduce
	 * @brief Base abstract class for reduce
	 */
	template <class LocalReduceType>
	class SimpleDynamicsReduce : public ParticleDynamicsReduce<LocalReduceType, size_t, ReduceRangeFunctor>
	{
	public:
		template <typename... Args>
		explicit SimpleDynamicsReduce(SPHBody &sph_body, Args &&...args)
			: ParticleDynamicsReduce<LocalReduceType, size_t, ReduceRangeFunctor>(
				  sph_body.BodyRange(), std::forward<Args>(args)...)
		{
			this->functor_ = std::bind(&LocalReduceType::reduceRange, this, _1, _2);
		};
		virtual ~SimpleDynamicsReduce(){};
	};

	/**
	 * @class OldParticleDynamicsReduce
	 * @brief Base abstract class for reduce
	 */
	template <class ReturnType, typename ReduceOperation>
	class OldParticleDynamicsReduce : public OldParticleDynamics<ReturnType>
	{
	public:
		explicit OldParticleDynamicsReduce(SPHBody &sph_body)
			: OldParticleDynamics<ReturnType>(sph_body), quantity_name_("ReducedQuantity"), initial_reference_(),
			  functor_reduce_function_(std::bind(&OldParticleDynamicsReduce::ReduceFunction, this, _1, _2)){};
		virtual ~OldParticleDynamicsReduce(){};

		ReturnType InitialReference() { return initial_reference_; };
		std::string QuantityName() { return quantity_name_; };

		virtual ReturnType exec(Real dt = 0.0) override
		{
			size_t all_real_particles = this->base_particles_->all_real_particles_;
			this->setBodyUpdated();
			SetupReduce();
			ReturnType temp = particle_reduce(all_real_particles,
											  initial_reference_, functor_reduce_function_, operation_, dt);
			return OutputResult(temp);
		};
		virtual ReturnType parallel_exec(Real dt = 0.0) override
		{
			size_t all_real_particles = this->base_particles_->all_real_particles_;
			this->setBodyUpdated();
			SetupReduce();
			ReturnType temp = particle_parallel_reduce(all_real_particles,
													   initial_reference_, functor_reduce_function_, operation_, dt);
			return this->OutputResult(temp);
		};

	protected:
		ReduceOperation operation_;
		std::string quantity_name_;

		/** initial or reference value */
		ReturnType initial_reference_;
		virtual void SetupReduce(){};
		virtual ReturnType ReduceFunction(size_t index_i, Real dt = 0.0) = 0;
		virtual ReturnType OutputResult(ReturnType reduced_value) { return reduced_value; };
		ReduceFunctor<ReturnType> functor_reduce_function_;
	};

	/**
	 * @class OldParticleDynamicsSimple
	 * @brief Simple particle dynamics without considering particle interaction
	 */
	class OldParticleDynamicsSimple : public OldParticleDynamics<void>
	{
	public:
		explicit OldParticleDynamicsSimple(SPHBody &sph_body)
			: OldParticleDynamics<void>(sph_body),
			  functor_update_(std::bind(&OldParticleDynamicsSimple::Update, this, _1, _2)){};
		virtual ~OldParticleDynamicsSimple(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		virtual void Update(size_t index_i, Real dt = 0.0) = 0;
		ParticleFunctor functor_update_;
	};

	/**
	 * @class OldInteractionDynamics
	 * @brief This is the class for particle interaction with other particles
	 */
	class OldInteractionDynamics : public OldParticleDynamics<void>
	{
	public:
		explicit OldInteractionDynamics(SPHBody &sph_body)
			: OldParticleDynamics<void>(sph_body),
			  functor_interaction_(std::bind(&OldInteractionDynamics::Interaction,
											 this, _1, _2)){};
		virtual ~OldInteractionDynamics(){};

		/** pre process such as update ghost state */
		StdVec<OldParticleDynamics<void> *> pre_processes_;
		/** post process such as impose constraint */
		StdVec<OldParticleDynamics<void> *> post_processes_;

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		friend class CombinedInteractionDynamics;
		virtual void Interaction(size_t index_i, Real dt = 0.0) = 0;
		ParticleFunctor functor_interaction_;
	};

	/**
	 * @class CombinedInteractionDynamics
	 * @brief This is the class for combining several interactions dynamics,
	 * which share the particle loop but are independent from each other,
	 * aiming to increase computing intensity under the data caching environment
	 */
	class CombinedInteractionDynamics : public OldInteractionDynamics
	{
	public:
		explicit CombinedInteractionDynamics(OldInteractionDynamics &dynamics_a, OldInteractionDynamics &dynamics_b);
		virtual ~CombinedInteractionDynamics(){};

	protected:
		OldInteractionDynamics &dynamics_a_, &dynamics_b_;
		virtual void setupDynamics(Real dt = 0.0) override;
		virtual void Interaction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	 * @class OldInteractionDynamicsWithUpdate
	 * @brief This class includes an interaction and a update steps
	 */
	class OldInteractionDynamicsWithUpdate : public OldInteractionDynamics
	{
	public:
		explicit OldInteractionDynamicsWithUpdate(SPHBody &sph_body)
			: OldInteractionDynamics(sph_body),
			  functor_update_(std::bind(&OldInteractionDynamicsWithUpdate::Update,
										this, _1, _2)) {}
		virtual ~OldInteractionDynamicsWithUpdate(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		virtual void Update(size_t index_i, Real dt = 0.0) = 0;
		ParticleFunctor functor_update_;
	};

	/**
	 * @class OldParticleDynamics1Level
	 * @brief This class includes an initialization, an interaction and a update steps
	 */
	class OldParticleDynamics1Level : public OldInteractionDynamicsWithUpdate
	{
	public:
		explicit OldParticleDynamics1Level(SPHBody &sph_body)
			: OldInteractionDynamicsWithUpdate(sph_body),
			  functor_initialization_(std::bind(&OldParticleDynamics1Level::Initialization,
												this, _1, _2)) {}
		virtual ~OldParticleDynamics1Level(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		virtual void Initialization(size_t index_i, Real dt = 0.0) = 0;
		ParticleFunctor functor_initialization_;
	};

	/**
	 * @class InteractionDynamicsSplitting
	 * @brief This is for the splitting algorithm
	 */
	class InteractionDynamicsSplitting : public OldInteractionDynamics
	{
	public:
		explicit InteractionDynamicsSplitting(SPHBody &sph_body);
		virtual ~InteractionDynamicsSplitting(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		SplitCellLists &split_cell_lists_;
	};
}
#endif // PARTICLE_DYNAMICS_ALGORITHMS_H