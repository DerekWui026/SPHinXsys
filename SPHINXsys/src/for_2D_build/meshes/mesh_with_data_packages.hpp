/**
 * @file 	base_mesh.hpp
 * @brief 	This is the implementation of the template function and class for base mesh
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef MESH_WITH_DATA_PACKAGES_2D_HPP
#define MESH_WITH_DATA_PACKAGES_2D_HPP

#include "mesh_with_data_packages.h"

//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	template <int PKG_SIZE, int ADDRS_SIZE>
	template <class DataType>
	DataType GridDataPackage<PKG_SIZE, ADDRS_SIZE>::
		probeDataPackage(PackageDataAddress<DataType> &pkg_data_addrs, const Vecd &position)
	{
		Vecu grid_idx = CellIndexFromPosition(position);
		Vecd grid_pos = GridPositionFromIndex(grid_idx);
		Vecd alpha = (position - grid_pos) / grid_spacing_;
		Vecd beta = Vec2d(1.0) - alpha;

		DataType bilinear = *pkg_data_addrs[grid_idx[0]][grid_idx[1]] * beta[0] * beta[1] +
							*pkg_data_addrs[grid_idx[0] + 1][grid_idx[1]] * alpha[0] * beta[1] +
							*pkg_data_addrs[grid_idx[0]][grid_idx[1] + 1] * beta[0] * alpha[1] +
							*pkg_data_addrs[grid_idx[0] + 1][grid_idx[1] + 1] * alpha[0] * alpha[1];

		return bilinear;
	}
	//=================================================================================================//
	template <int PKG_SIZE, int ADDRS_SIZE>
	template <typename InDataType, typename OutDataType>
	void GridDataPackage<PKG_SIZE, ADDRS_SIZE>::
		computeGradient(PackageDataAddress<InDataType> &in_pkg_data_addrs,
						PackageDataAddress<OutDataType> out_pkg_data_addrs, Real dt)
	{
		for (int i = 1; i != PKG_SIZE + 1; ++i)
			for (int j = 1; j != PKG_SIZE + 1; ++j)
			{
				Real dphidx = (*in_pkg_data_addrs[i + 1][j] - *in_pkg_data_addrs[i - 1][j]);
				Real dphidy = (*in_pkg_data_addrs[i][j + 1] - *in_pkg_data_addrs[i][j - 1]);
				*out_pkg_data_addrs[i][j] = 0.5 * Vecd(dphidx, dphidy) / grid_spacing_;
			}
	}
	//=================================================================================================//
	template <int PKG_SIZE, int ADDRS_SIZE>
	template <typename InDataType, typename OutDataType>
	void GridDataPackage<PKG_SIZE, ADDRS_SIZE>::
		computeNormalizedGradient(PackageDataAddress<InDataType> &in_pkg_data_addrs,
								  PackageDataAddress<OutDataType> out_pkg_data_addrs, Real dt)
	{
		for (int i = 1; i != PKG_SIZE + 1; ++i)
			for (int j = 1; j != PKG_SIZE + 1; ++j)
			{
				Real dphidx = (*in_pkg_data_addrs[i + 1][j] - *in_pkg_data_addrs[i - 1][j]);
				Real dphidy = (*in_pkg_data_addrs[i][j + 1] - *in_pkg_data_addrs[i][j - 1]);
				Vecd normal = Vecd(dphidx, dphidy);
				*out_pkg_data_addrs[i][j] = normal / (normal.norm() + TinyReal);
			}
	}
	//=================================================================================================//
	template <int PKG_SIZE, int ADDRS_SIZE>
	template <typename DataType>
	void GridDataPackage<PKG_SIZE, ADDRS_SIZE>::initializePackageDataAddress<DataType>::
	operator()(GeneralDataAssemble<PackageData> &all_pkg_data,
			   GeneralDataAssemble<PackageDataAddress> &all_pkg_data_addrs)
	{
		constexpr int type_index = DataTypeIndex<DataType>::value;
		for (size_t l = 0; l != std::get<type_index>(all_pkg_data).size(); ++l)
		{
			PackageData<DataType> &pkg_data = *std::get<type_index>(all_pkg_data)[l];
			PackageDataAddress<DataType> &pkg_data_addrs = *std::get<type_index>(all_pkg_data_addrs)[l];
			for (int i = 0; i != ADDRS_SIZE; ++i)
				for (int j = 0; j != ADDRS_SIZE; ++j)
				{
					pkg_data_addrs[i][j] = &pkg_data[0][0];
				}
		}
	}
	//=================================================================================================//
	template <int PKG_SIZE, int ADDRS_SIZE>
	template <typename DataType>
	DataType GridDataPackage<PKG_SIZE, ADDRS_SIZE>::
		CornerAverage(PackageDataAddress<DataType> &pkg_data_addrs, Veci addrs_index, Veci corner_direction)
	{
		DataType average(0);
		for (int i = 0; i != 2; ++i)
			for (int j = 0; j != 2; ++j)
			{
				int x_index = addrs_index[0] + i * corner_direction[0];
				int y_index = addrs_index[1] + j * corner_direction[1];
				average += *pkg_data_addrs[x_index][y_index];
			}
		return average * 0.25;
	}
	//=================================================================================================//
	template <int PKG_SIZE, int ADDRS_SIZE>
	template <typename DataType>
	void GridDataPackage<PKG_SIZE, ADDRS_SIZE>::assignPackageDataAddress<DataType>::
	operator()(GeneralDataAssemble<PackageDataAddress> &all_pkg_data_addrs,
			   const Vecu &addrs_index,
			   GeneralDataAssemble<PackageData> &all_pkg_data,
			   const Vecu &data_index)
	{
		constexpr int type_index = DataTypeIndex<DataType>::value;
		for (size_t l = 0; l != std::get<type_index>(all_pkg_data).size(); ++l)
		{
			PackageData<DataType> &pkg_data = *std::get<type_index>(all_pkg_data)[l];
			PackageDataAddress<DataType> &pkg_data_addrs = *std::get<type_index>(all_pkg_data_addrs)[l];
			pkg_data_addrs[addrs_index[0]][addrs_index[1]] = &pkg_data[data_index[0]][data_index[1]];
		}
	} //=================================================================================================//
	template <class MeshFieldType, class GridDataPackageType>
	template <typename DataType, typename PackageDataType, PackageDataType GridDataPackageType::*MemPtr>
	DataType MeshWithGridDataPackages<MeshFieldType, GridDataPackageType>::
		DataValueFromGlobalIndex(const Vecu &global_grid_index)
	{
		Vecu pkg_index_(0);
		Vecu local_data_index(0);
		for (int n = 0; n != 2; n++)
		{
			size_t cell_index_in_this_direction = global_grid_index[n] / pkg_size_;
			pkg_index_[n] = cell_index_in_this_direction;
			local_data_index[n] = global_grid_index[n] - cell_index_in_this_direction * pkg_size_;
		}
		PackageDataType &data = data_pkg_addrs_[pkg_index_[0]][pkg_index_[1]]->*MemPtr;
		return data[local_data_index[0]][local_data_index[1]];
	}
	//=================================================================================================//
	template <class MeshFieldType, class GridDataPackageType>
	void MeshWithGridDataPackages<MeshFieldType, GridDataPackageType>::
		initializePackageAddressesInACell(const Vecu &cell_index)
	{
		int i = (int)cell_index[0];
		int j = (int)cell_index[1];

		GridDataPackageType *data_pkg = data_pkg_addrs_[i][j];
		if (data_pkg->is_inner_pkg_)
		{
			for (int l = 0; l != pkg_addrs_size_; ++l)
				for (int m = 0; m != pkg_addrs_size_; ++m)
				{
					std::pair<int, int> x_pair = CellShiftAndDataIndex(l);
					std::pair<int, int> y_pair = CellShiftAndDataIndex(m);
					data_pkg->assignAllPackageDataAddress(
						Vecu(l, m),
						data_pkg_addrs_[i + x_pair.first][j + y_pair.first],
						Vecu(x_pair.second, y_pair.second));
				}
		}
	}
	//=================================================================================================//
	template <class MeshFieldType, class GridDataPackageType>
	void MeshWithGridDataPackages<MeshFieldType, GridDataPackageType>::allocateMeshDataMatrix()
	{
		Allocate2dArray(data_pkg_addrs_, number_of_cells_);
	}
	//=================================================================================================//
	template <class MeshFieldType, class GridDataPackageType>
	void MeshWithGridDataPackages<MeshFieldType, GridDataPackageType>::deleteMeshDataMatrix()
	{
		Delete2dArray(data_pkg_addrs_, number_of_cells_);
	}
	//=================================================================================================//
	template <class MeshFieldType, class GridDataPackageType>
	void MeshWithGridDataPackages<MeshFieldType, GridDataPackageType>::
		assignDataPackageAddress(const Vecu &cell_index, GridDataPackageType *data_pkg)
	{
		data_pkg_addrs_[cell_index[0]][cell_index[1]] = data_pkg;
	}
	//=================================================================================================//
	template <class MeshFieldType, class GridDataPackageType>
	GridDataPackageType *MeshWithGridDataPackages<MeshFieldType, GridDataPackageType>::
		DataPackageFromCellIndex(const Vecu &cell_index)
	{
		return data_pkg_addrs_[cell_index[0]][cell_index[1]];
	}
	//=================================================================================================//
	template <class MeshFieldType, class GridDataPackageType>
	template <class DataType, typename PackageDataAddressType, PackageDataAddressType GridDataPackageType::*MemPtr>
	DataType MeshWithGridDataPackages<MeshFieldType, GridDataPackageType>::probeMesh(const Vecd &position)
	{
		Vecu grid_index = CellIndexFromPosition(position);
		size_t i = grid_index[0];
		size_t j = grid_index[1];

		GridDataPackageType *data_pkg = data_pkg_addrs_[i][j];
		PackageDataAddressType &pkg_data_addrs = data_pkg->*MemPtr;
		return data_pkg->is_inner_pkg_ ? data_pkg->GridDataPackageType::template probeDataPackage<DataType>(pkg_data_addrs, position)
									   : *pkg_data_addrs[0][0];
	}
	//=================================================================================================//
}
//=================================================================================================//
#endif // MESH_WITH_DATA_PACKAGES_2D_HPP