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
 * @file 	base_particles.h
 * @brief 	This is the base class of SPH particles. The basic data of the particles
 *			is saved in separated large vectors. Each derived class will introduce several extra
 * 			vectors for the new data. Note that there is no class of single particle.
 * @author	Xiangyu Hu and Chi Zhang
 */

#ifndef BASE_PARTICLES_H
#define BASE_PARTICLES_H

#include "base_data_package.h"
#include "sph_data_containers.h"
#include "base_material.h"
#include "xml_engine.h"

#include <fstream>

namespace SPH
{

	class SPHBody;
	class BaseMaterial;
	class ParticleGenerator;
	class BodySurface;
	template <class ReturnType>
	class BaseParticleDynamics;

	/**
	 * @class BaseParticles
	 * @brief Particles with essential (geometric and kinematic) data.
	 * There are three types of particles， all particles of a same type are saved with continuous memory segments.
	 * The first type is real particles whose states are updated by particle dynamics.
	 * One is buffer particles whose state are not updated by particle dynamics.
	 * Buffer particles are saved behind real particles.
	 * The global value of total_real_particles_ separate the real and buffer particles.
	 * They may be switched from real particles or switch to real particles.
	 * As the memory for both particles are continuous, such switch is achieved at the memory boundary sequentially.
	 * The basic idea is swap the data of the last real particle with the one will be switched particle,
	 * and then switch this swapped last particle as buffer particle by decrease the total_real_particles_ by one.
	 * Switch from buffer particle to real particle is easy. One just need to assign expect state to
	 * the first buffer particle and increase total_real_particles_ by one.
	 * The other is ghost particles whose states are updated according to
	 * boundary condition if their indices are included in the neighbor particle list.
	 * The ghost particles are saved behind the buffer particles.
	 * The global value of real_particles_bound_ separate the sum of real and buffer particles with ghost particles.
	 * The global value of total_ghost_particles_ indicates the total number of ghost particles in use.
	 * It will be initialized to zero before a time step.
	 * In SPHinXsys, the variables registered in general particle data (ParticleData) belong to a heirachy of two layers.
	 * The first is for the global basic physical states to describe the physical process.
	 * These variables are defined within the classes of particles.
	 * The second is for the local, dynamics-method-related variables, which are defined in specific methods,
	 * and are only used by the relevant methods.
	 * There is a rule of single registration, that is,
	 * a variable is only allowed to be registered with a name once by the function registerAVariable.
	 * The usage of the second- and third-layer variables is accessed by getVariableByName.
	 * Such a rule requires careful design of the code.
	 */
	class BaseParticles
	{
	private:
		UniquePtrKeepers<BaseParticleDynamics<void>> derived_particle_data_;

	public:
		explicit BaseParticles(SPHBody &sph_body, BaseMaterial *base_material);
		virtual ~BaseParticles(){};

		StdLargeVec<Vecd> pos_n_;		  /**< current position */
		StdLargeVec<Vecd> vel_n_;		  /**< current particle velocity */
		StdLargeVec<Vecd> dvel_dt_;		  /**< total acceleration including inner pressure- or stress-induced acceleration and other accelerations */
		StdLargeVec<Vecd> dvel_dt_prior_; /**< other, such as gravity and viscous, accelerations */

		StdLargeVec<Real> Vol_;	  /**< particle volume */
		StdLargeVec<Real> rho_n_; /**< current particle density */
		StdLargeVec<Real> mass_;  /**< particle mass */
		//----------------------------------------------------------------------
		// Global information for all particles
		//----------------------------------------------------------------------
		Real rho0_;				/**< reference density*/
		Real sigma0_;			/**< reference number density. */
		Real speed_max_;		/**< Maximum particle speed. */
		Real signal_speed_max_; /**< Maximum signal speed.*/
		//----------------------------------------------------------------------
		// Global information for defining particle groups
		//----------------------------------------------------------------------
		size_t total_real_particles_;
		size_t real_particles_bound_; /**< Maximum possible number of real particles. Also the start index of ghost particles. */
		size_t total_ghost_particles_;
		//----------------------------------------------------------------------
		//		Generalized particle data for parameterized management
		//----------------------------------------------------------------------
		ParticleData all_particle_data_;
		ParticleDataMap all_variable_maps_;
		StdVec<BaseParticleDynamics<void> *> derived_variables_;
		ParticleVariableList variables_to_write_;

		/** register a variable defined in a class (can be non-particle class) */
		template <typename VariableType>
		void registerAVariable(StdLargeVec<VariableType> &variable_addrs,
							   const std::string &variable_name, VariableType initial_value = VariableType(0));

		/** register a variable and copying data from an exist variable */
		template <typename VariableType>
		void registerAVariable(StdLargeVec<VariableType> &variable_addrs,
							   const std::string &new_variable_name, const std::string &old_variable_name);

		/** get a registered variable from particles by its name. return by pointer so that return nullptr if fail. */
		template <typename VariableType>
		StdLargeVec<VariableType> *getVariableByName(const std::string &variable_name);

		/** add a variable into a particle vairable name list */
		template <typename VariableType>
		void addAVariableNameToList(ParticleVariableList &variable_name_list, const std::string &variable_name);

		/** add a variable into the list for state output */
		template <typename VariableType>
		void addAVariableToWrite(const std::string &variable_name);

		/** add a derived variable into the list for state output */
		template <class DerivedVariableMethod>
		void addDerivedVariableToWrite();

		/** add a variable into the list for restart */
		template <typename VariableType>
		void addAVariableToRestart(const std::string &variable_name);

		/** add a variable into the list for particle reload */
		template <typename VariableType>
		void addAVariableToReload(const std::string &variable_name);

		//----------------------------------------------------------------------
		//		Particle data for sorting
		//----------------------------------------------------------------------
		StdLargeVec<size_t> unsorted_id_; /**< the ids assigned just after particle generated. */
		StdLargeVec<size_t> sorted_id_;	/**< the sorted particle ids of particles from unsorted ids. */
		StdLargeVec<size_t> sequence_;	/**< the sequence referred for sorting. */
		ParticleData sortable_data_;
		ParticleDataMap sortable_variable_maps_;

		/** register an already defined variable as sortable */
		template <typename VariableType>
		void registerASortableVariable(const std::string &variable_name);

		SPHBody *getSPHBody() { return sph_body_; };
		/** initialize other variables  based one geometric variables and material */
		virtual void initializeOtherVariables();
		void addBufferParticles(size_t buffer_size);
		void copyFromAnotherParticle(size_t this_index, size_t another_index);
		void updateFromAnotherParticle(size_t this_index, size_t another_index);
		size_t insertAGhostParticle(size_t index_i);
		void switchToBufferParticle(size_t index_i);

		/** Write particle data in Vtk format for Paraview. */
		template <typename OutStreamType>
		void writeParticlesToVtk(OutStreamType &output_stream);
		/** Write particle data in PLT format for Tecplot. */
		void writeParticlesToPltFile(std::ofstream &output_file);
		/** Write only surface particle data in VTU format for Paraview. TODO: this should be generalized for body part by particles */
		virtual void writeSurfaceParticlesToVtuFile(std::ostream &output_file, BodySurface &surface_particles);

		void resizeXmlDocForParticles(XmlEngine &xml_engine);
		void writeParticlesToXmlForRestart(std::string &filefullpath);
		void readParticleFromXmlForRestart(std::string &filefullpath);
		XmlEngine *getReloadXmlEngine() { return &reload_xml_engine_; };
		void writeToXmlForReloadParticle(std::string &filefullpath);
		void readFromXmlForReloadParticle(std::string &filefullpath);

		virtual BaseParticles *ThisObjectPtr() { return this; };

		/** Normalize the kernel gradient. */
		virtual Vecd normalizeKernelGradient(size_t particle_index_i, Vecd &kernel_gradient) { return kernel_gradient; };
		/** Get the kernel gradient in weak form. */
		virtual Vecd getKernelGradient(size_t particle_index_i, size_t particle_index_j,
									   Real dW_ij, Vecd &e_ij)
		{
			return dW_ij * e_ij;
		};

	protected:
		SPHBody *sph_body_; /**< The body in which the particles belongs to. */
		std::string body_name_;
		XmlEngine restart_xml_engine_;
		XmlEngine reload_xml_engine_;
		ParticleVariableList variables_to_restart_;
		ParticleVariableList variables_to_reload_;
		void addAParticleEntry();

		virtual void writePltFileHeader(std::ofstream &output_file);
		virtual void writePltFileParticleData(std::ofstream &output_file, size_t index_i);

		/** resize a particle data. */
		template <typename VariableType>
		struct resizeParticleData
		{
			void operator()(ParticleData &particle_data, size_t new_size) const;
		};

		/** Fill a particle variable with default data. */
		template <typename VariableType>
		struct addAParticleDataValue
		{
			void operator()(ParticleData &particle_data) const;
		};

		/** Copy a particle variable value from another particle. */
		template <typename VariableType>
		struct copyAParticleDataValue
		{
			void operator()(ParticleData &particle_data, size_t this_index, size_t another_index) const;
		};

		ParticleDataOperation<resizeParticleData> resize_particle_data_;
		ParticleDataOperation<addAParticleDataValue> add_a_particle_value_;
		ParticleDataOperation<copyAParticleDataValue> copy_a_particle_value_;
	};

	struct WriteAParticleVariableToXml
	{
		XmlEngine &xml_engine_;
		size_t &total_real_particles_;
		WriteAParticleVariableToXml(XmlEngine &xml_engine, size_t &total_real_particles)
			: xml_engine_(xml_engine), total_real_particles_(total_real_particles){};

		template <typename VariableType>
		void operator()(std::string &variable_name, StdLargeVec<VariableType> &variable) const;
	};

	struct ReadAParticleVariableFromXml
	{
		XmlEngine &xml_engine_;
		size_t &total_real_particles_;
		ReadAParticleVariableFromXml(XmlEngine &xml_engine, size_t &total_real_particles)
			: xml_engine_(xml_engine), total_real_particles_(total_real_particles){};

		template <typename VariableType>
		void operator()(std::string &variable_name, StdLargeVec<VariableType> &variable) const;
	};

	/**
	 * @class BaseDerivedVariable
	 * @brief computing displacement from current and initial particle position
	 */
	template <typename VariableType>
	class BaseDerivedVariable
	{
	public:
		using DerivedVariableType = VariableType;
		std::string variable_name_;

		BaseDerivedVariable(const SPHBody &sph_body, const std::string &variable_name);
		virtual ~BaseDerivedVariable(){};

	protected:
		StdLargeVec<VariableType> derived_variable_;
	};
}
#endif // BASE_PARTICLES_H
