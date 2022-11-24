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
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
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
 * @file 	mesh_with_data_packages.h
 * @brief 	This class is designed to save memory and increase computational efficiency
 *	on mesh. //TODO: the connection between successive meshes in refined mesh should enhanced.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef MESH_WITH_DATA_PACKAGES_H
#define MESH_WITH_DATA_PACKAGES_H

#include "base_mesh.h"
#include "base_variables.h"
#include "my_memory_pool.h"

#include <fstream>
#include <algorithm>
#include <mutex>
#include <functional>
using namespace std::placeholders;

namespace SPH
{
	/** Iterator on a collection of mesh data packages. sequential computing. */
	template <class DataPackageType, typename LocalFunction, typename... Args>
	void package_for(const ConcurrentVec<DataPackageType *> &data_pkgs,
					 const LocalFunction &local_function, Args &&...args)
	{
		for (size_t i = 0; i != data_pkgs.size(); ++i)
			local_function(i);
	};
	/** Iterator on a collection of mesh data packages. parallel computing. */
	template <class DataPackageType, typename LocalFunction, typename... Args>
	void package_parallel_for(const ConcurrentVec<DataPackageType *> &data_pkgs,
							  const LocalFunction &local_function, Args &&...args)
	{
		parallel_for(
			blocked_range<size_t>(0, data_pkgs.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					local_function(i);
				}
			},
			ap);
	};

	/**
	 * @class BaseDataPackage
	 * @brief Abstract base class for a data package,
	 * by which the data in a derived class can be on- or off-grid.
	 * The data package can be defined in a cell of a background mesh so the pkg_index is
	 * the cell location on the mesh.
	 * TODO: The class will be enriched with general methods for all data packages.
	 */
	class BaseDataPackage
	{
	public:
		Vecu pkg_index_;	/**< index of this data package on the background mesh, Vecu(0) if it is not on the mesh. */
		bool is_inner_pkg_; /**< If true, its data package is on the background mesh. */
		/** reserved value: 0 not occupying background mesh, 1 occupying.
		 *  guide to use: large magnitude for high priority of the data package.
		 */
		int state_indicator_;

		BaseDataPackage() : pkg_index_(0), is_inner_pkg_(false), state_indicator_(0){};
		virtual ~BaseDataPackage(){};
	};

	/**
	 * @class GridDataPackage
	 * @brief Abstract base class for a data package
	 * whose data are defined on the grids of a small mesh patch.
	 * note tha pkg_addrs_size = pkg_size + 2 * pkg_addrs_buffer;
	 * Also note that, while the mesh lower bound locates the first data address,
	 * the data lower bound locates the first data.
	 * The package lower bound is just in the middle of mesh and data lower bounds.
	 */
	template <int PKG_SIZE, int ADDRS_BUFFER>
	class GridDataPackage : public BaseDataPackage, public BaseMesh
	{
	public:
		static constexpr int pkg_size = PKG_SIZE;
		static constexpr int pkg_addrs_size = PKG_SIZE + ADDRS_BUFFER * 2;
		static constexpr int pkg_addrs_buffer = ADDRS_BUFFER;
		static constexpr int pkg_ops_end = PKG_SIZE + pkg_addrs_buffer;
		template <typename DataType>
		using PackageData = PackageDataMatrix<DataType, PKG_SIZE>;
		template <typename DataType>
		using PackageDataAddress = PackageDataMatrix<DataType *, pkg_addrs_size>;
		/** Matrix data for temporary usage. Note that it is array with pkg_addrs_size.  */
		template <typename DataType>
		using PackageTemporaryData = PackageDataMatrix<DataType, pkg_addrs_size>;

		DataContainerAddressAssemble<PackageData> all_pkg_data_;
		DataContainerAddressAssemble<PackageDataAddress> all_pkg_data_addrs_;
		DataContainerAssemble<PackageData> extra_pkg_data_;
		DataContainerAssemble<PackageDataAddress> extra_pkg_data_addrs_;

		GridDataPackage() : BaseDataPackage(), BaseMesh(Vecu(pkg_addrs_size)){};
		virtual ~GridDataPackage(){};
		/** lower bound coordinate for the data as reference */
		Vecd DataLowerBound() { return mesh_lower_bound_ + Vecd(grid_spacing_) * (Real)pkg_addrs_buffer; };
		/** initialize package mesh geometric information. */
		void initializePackageGeometry(const Vecd &pkg_lower_bound, Real data_spacing)
		{
			mesh_lower_bound_ = pkg_lower_bound - Vecd(data_spacing) * ((Real)pkg_addrs_buffer - 0.5);
			grid_spacing_ = data_spacing;
		};
		/** void (non_value_returning) function iterate on all data points by value,
		 *  for function only involving the data itself. */
		template <typename FunctionOnData>
		void for_each_data(const FunctionOnData &function);
		/** void (non_value_returning) function iterate on all data points by address,
		 *  for function involving operations on data neighbors. */
		template <typename FunctionOnAddress>
		void for_each_addrs(const FunctionOnAddress &function);

		// upwind algorithm choosing candidate difference by the sign
		Real upwindDifference(Real sign, Real df_p, Real df_n)
		{
			if (sign * df_p >= 0.0 && sign * df_n >= 0.0)
				return df_n;
			if (sign * df_p <= 0.0 && sign * df_n <= 0.0)
				return df_p;
			if (sign * df_p > 0.0 && sign * df_n < 0.0)
				return 0.0;

			Real df = df_p;
			if (sign * df_p < 0.0 && sign * df_n > 0.0)
			{
				Real ss = sign * (fabs(df_p) - fabs(df_n)) / (df_p - df_n);
				if (ss > 0.0)
					df = df_n;
			}

			return df;
		};

		/** access specific package data with discrete variable */
		template <typename DataType>
		PackageData<DataType> &getPackageData(const DiscreteVariable<DataType> &discrete_variable)
		{
			constexpr int type_index = DataTypeIndex<DataType>::value;
			return std::get<type_index>(extra_pkg_data_)[discrete_variable.IndexInContainer()];
		};
		/** access specific package data with discrete variable */
		template <typename DataType>
		PackageDataAddress<DataType> &getPackageDataAddress(const DiscreteVariable<DataType> &discrete_variable)
		{
			constexpr int type_index = DataTypeIndex<DataType>::value;
			return std::get<type_index>(extra_pkg_data_addrs_)[discrete_variable.IndexInContainer()];
		};
		/** probe by applying bi and tri-linear interpolation within the package. */
		template <typename DataType>
		DataType probeDataPackage(PackageDataAddress<DataType> &pkg_data_addrs, const Vecd &position);
		/** compute gradient transform within data package */
		template <typename InDataType, typename OutDataType>
		void computeGradient(PackageDataAddress<InDataType> &in_pkg_data_addrs,
							 PackageDataAddress<OutDataType> out_pkg_data_addrs, Real dt = 0.0);
		/** assign value to data package according to grid position */
		template <typename DataType, typename FunctionByPosition>
		void assignByPosition(const DiscreteVariable<DataType> &discrete_variable,
							  const FunctionByPosition &function_by_position);
		/** compute gradient transform within data package */
		template <typename InDataType, typename OutDataType>
		void computeGradient(const DiscreteVariable<InDataType> &in_variable, const DiscreteVariable<OutDataType> &out_variable);

	protected:
		/** register a variable defined in a class (can be non-particle class) */
		template <typename DataType>
		void registerPackageData(PackageData<DataType> &pkg_data,
								 PackageDataAddress<DataType> &pkg_data_addrs)
		{
			constexpr int type_index = DataTypeIndex<DataType>::value;
			std::get<type_index>(all_pkg_data_).push_back(&pkg_data);
			std::get<type_index>(all_pkg_data_addrs_).push_back(&pkg_data_addrs);
		};

		/** set the initial package data address within a derived class constructor */
		template <typename DataType>
		struct initializePackageDataAddress
		{
			void operator()(DataContainerAddressAssemble<PackageData> &all_pkg_data,
							DataContainerAddressAssemble<PackageDataAddress> &all_pkg_data_addrs);
		};
		DataAssembleOperation<initializePackageDataAddress> initialize_pkg_data_addrs_;

		/** set the initial package data address within a derived class constructor */
		template <typename DataType>
		struct initializeExtraPackageDataAddress
		{
			void operator()(DataContainerAssemble<PackageData> &all_pkg_data,
							DataContainerAssemble<PackageDataAddress> &all_pkg_data_addrs);
		};
		DataAssembleOperation<initializeExtraPackageDataAddress> initialize_extra_pkg_data_addrs_;

		/** assign address for a package data when the package is an inner one */
		template <typename DataType>
		struct assignPackageDataAddress
		{
			void operator()(DataContainerAddressAssemble<PackageDataAddress> &all_pkg_data_addrs,
							const Vecu &addrs_index,
							DataContainerAddressAssemble<PackageData> &all_pkg_data,
							const Vecu &data_index);
		};
		DataAssembleOperation<assignPackageDataAddress> assign_pkg_data_addrs_;

		/** assign address for extra package data when the package is an inner one */
		template <typename DataType>
		struct assignExtraPackageDataAddress
		{
			void operator()(DataContainerAssemble<PackageDataAddress> &all_pkg_data_addrs,
							const Vecu &addrs_index,
							DataContainerAssemble<PackageData> &all_pkg_data,
							const Vecu &data_index);
		};
		DataAssembleOperation<assignExtraPackageDataAddress> assign_extra_pkg_data_addrs_;

		/** allocate memory for extra package data when the package is an inner one */
		template <typename DataType>
		struct ExtraVariablesAllocation
		{
			void operator()(DataContainerAssemble<PackageData> &extra_pkg_data,
							DataContainerAssemble<PackageDataAddress> &extra_pkg_data_addrs,
							const DiscreteVariableAssemble &extra_variables)
			{
				constexpr int type_index = DataTypeIndex<DataType>::value;
				size_t total_variables = std::get<type_index>(extra_variables).size();
				std::get<type_index>(extra_pkg_data).resize(total_variables);
				std::get<type_index>(extra_pkg_data_addrs).resize(total_variables);
			};
		};
		DataAssembleOperation<ExtraVariablesAllocation> allocate_extra_variables_;

		/** obtain averaged value at a corner of a data cell */
		template <typename DataType>
		DataType CornerAverage(PackageDataAddress<DataType> &pkg_data_addrs, Veci addrs_index, Veci corner_direction);

	public:
		void initializeSingularDataAddress()
		{
			initialize_pkg_data_addrs_(all_pkg_data_, all_pkg_data_addrs_);
			initialize_extra_pkg_data_addrs_(extra_pkg_data_, extra_pkg_data_addrs_);
		};

		void assignAllPackageDataAddress(const Vecu &addrs_index, GridDataPackage *src_pkg, const Vecu &data_index)
		{
			assign_pkg_data_addrs_(all_pkg_data_addrs_, addrs_index, src_pkg->all_pkg_data_, data_index);
			assign_extra_pkg_data_addrs_(extra_pkg_data_addrs_, addrs_index, src_pkg->extra_pkg_data_, data_index);
		};

		void allocateExtraVariables(const DiscreteVariableAssemble &extra_variables)
		{
			allocate_extra_variables_(extra_pkg_data_, extra_pkg_data_addrs_, extra_variables);
		};
	};

	/**
	 * @class MeshWithGridDataPackages
	 * @brief Abstract class for mesh with grid-based data packages.
	 * @details The idea is to save sparse data on a cell-based mesh.
	 * We say sparse data, it means that only in some inner mesh cells there are no trivial data.
	 * A typical example is a level-set field which only has meaningful values near the interface,
	 * while the latter is in the inner region of a mesh.
	 * In this class, only some inner mesh cells are filled with data packages.
	 * Each data package is again a mesh, but grid based, where two sets of data are saved on its grid points.
	 * One is the field data of matrices with pkg_size, the other is corresponding address data of matrices with pkg_addrs_size.
	 * For two neighboring data packages, they share the data in the buffer which is in the overlap region.
	 * The filling of field data is achieved first by the data matrices by the function initializeDataInACell
	 * and then the address matrix by the function initializeAddressesInACell.
	 * All these data packages are indexed by a concurrent vector inner_data_pkgs_.
	 * Note that a data package should be not near the mesh bound, otherwise one will encounter the error "out of range".
	 */
	template <class MeshFieldType, class GridDataPackageType>
	class MeshWithGridDataPackages : public MeshFieldType, public Mesh
	{
	public:
		MyMemoryPool<GridDataPackageType> data_pkg_pool_;	   /**< memory pool for all packages in the mesh. */
		MeshDataMatrix<GridDataPackageType *> data_pkg_addrs_; /**< Address of data packages. */
		ConcurrentVec<GridDataPackageType *> inner_data_pkgs_; /**< Inner data packages which is able to carry out spatial operations. */

		template <typename... Args>
		explicit MeshWithGridDataPackages(BoundingBox tentative_bounds, Real data_spacing, size_t buffer_size, Args &&...args)
			: MeshFieldType(std::forward<Args>(args)...),
			  Mesh(tentative_bounds, GridDataPackageType::pkg_size * data_spacing, buffer_size),
			  data_spacing_(data_spacing),
			  global_mesh_(this->mesh_lower_bound_ + Vecd(data_spacing) * 0.5, data_spacing, this->number_of_cells_ * pkg_size)
		{
			allocateMeshDataMatrix();
		};
		virtual ~MeshWithGridDataPackages() { deleteMeshDataMatrix(); };

		void allocateMeshDataMatrix(); /**< allocate memories for addresses of data packages. */
		void deleteMeshDataMatrix();   /**< delete memories for addresses of data packages. */

		/** This function probe a mesh value */
		template <class DataType, typename PackageDataAddressType, PackageDataAddressType GridDataPackageType::*MemPtr>
		DataType probeMesh(const Vecd &position);
		template <class DataType>
		DataType probeMesh(const DiscreteVariable<DataType> &discrete_variable, const Vecd &position);
		/** spacing between the data, which is 1/ pkg_size of this grid spacing */
		virtual Real DataSpacing() override { return data_spacing_; };

	protected:
		const Real data_spacing_;													   /**< spacing of data in the data packages*/
		static constexpr int pkg_size = GridDataPackageType::pkg_size;				   /**< the size of the data package matrix*/
		static constexpr int pkg_addrs_buffer = GridDataPackageType::pkg_addrs_buffer; /**< the size of address buffer, a value less than the package size. */
		static constexpr int pkg_ops_end = GridDataPackageType::pkg_ops_end;		   /**< the size of operation loops. */
		static constexpr int pkg_addrs_size = GridDataPackageType::pkg_addrs_size;	   /**< the size of address matrix in the data packages. */
		std::mutex mutex_my_pool;													   /**< mutex exclusion for memory pool */
		BaseMesh global_mesh_;														   /**< the mesh for the locations of all possible data points. */
		/** Singular data packages. provided for far field condition with usually only two values.
		 * For example, when level set is considered. The first value for inner far-field and second for outer far-field */
		StdVec<GridDataPackageType *> singular_data_pkgs_addrs_;

		template <typename InitializeSingularData>
		void initializeASingularDataPackage(
			const DataContainerAddressAssemble<DiscreteVariable> &all_variables,
			const InitializeSingularData &initialize_singular_data)
		{
			GridDataPackageType *new_data_pkg = data_pkg_pool_.malloc();
			new_data_pkg->registerAllVariables();
			new_data_pkg->allocateExtraVariables(all_variables);
			initialize_singular_data(*new_data_pkg);
			new_data_pkg->initializeSingularDataAddress();
			singular_data_pkgs_addrs_.push_back(new_data_pkg);
		};

		void assignDataPackageAddress(const Vecu &cell_index, GridDataPackageType *data_pkg);
		GridDataPackageType *DataPackageFromCellIndex(const Vecu &cell_index);
		/** This function find the value of data from its index from global mesh. */
		template <typename DataType, typename PackageDataType, PackageDataType GridDataPackageType::*MemPtr>
		DataType DataValueFromGlobalIndex(const Vecu &global_grid_index);
		template <typename DataType>
		DataType DataValueFromGlobalIndex(const DiscreteVariable<DataType> &discrete_variable,
										  const Vecu &global_grid_index);

		void initializePackageAddressesInACell(const Vecu &cell_index);
		/** find related cell index and data index for a data package address matrix */
		std::pair<int, int> CellShiftAndDataIndex(int data_addrs_index_component)
		{
			std::pair<int, int> shift_and_index;
			int signed_date_index = data_addrs_index_component - pkg_addrs_buffer;
			shift_and_index.first = (signed_date_index + pkg_size) / pkg_size - 1;
			shift_and_index.second = signed_date_index - shift_and_index.first * pkg_size;
			return shift_and_index;
		}
	};
}
#endif // MESH_WITH_DATA_PACKAGES_H