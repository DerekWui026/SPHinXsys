/**
 * @file 	generative_structures.cpp
 * @author	Xiangyu Hu
 */

#include "complex_body.h"

#include "base_material.h"
#include "base_particles.h"
#include "neighbor_relation.h"
#include "adaptation.h"
#include "base_particle_dynamics.h"

namespace SPH
{
	//=================================================================================================//
	TreeBody::TreeBody(SPHSystem &sph_system, SharedPtr<Shape> shape_ptr)
		: SecondaryStructure(), RealBody(sph_system, shape_ptr), last_branch_id_(0)
	{
		root_ =  branches_ptr_keeper_.createPtr<Branch>(this);
	}
	//=================================================================================================//
	void TreeBody::buildParticleConfiguration(ParticleConfiguration &particle_configuration)
	{
		/** First branch
		 * Note that the first branch has only one particle.
		 * Find the neighbors in child branch, the first branch only have one child, id = 1.
		 */
		size_t particle_id = branches_[0]->inner_particles_.front();
		std::vector<size_t> neighboring_ids;
		neighboring_ids.push_back(branches_[1]->inner_particles_[0]);
		neighboring_ids.push_back(branches_[1]->inner_particles_[1]);
		/** Build configuration. */
		const StdLargeVec<Vecd> &pos_ =  base_particles_->pos_;
		NeighborRelationInner neighbor_relation_inner(this);
		for (size_t n = 0; n != neighboring_ids.size(); ++n)
		{
			Vecd displacement = pos_[particle_id] - pos_[neighboring_ids[n]];
			Neighborhood &neighborhood = particle_configuration[particle_id];
			neighbor_relation_inner(neighborhood, displacement, particle_id, neighboring_ids[n]);
		}
		/** Second branch. 
		 * The second branch has special parent branch, branch 0, consisting only one point.
		 * The child branch are two normal branch. 
		 */
		size_t num_ele = branches_[1]->inner_particles_.size();
		std::vector<size_t> child_ids;
		for (size_t k = 0; k < branches_[1]->out_edge_.size(); ++k)
		{
			child_ids.push_back(branches_[1]->out_edge_[k]);
		}

		for (size_t i = 0; i != num_ele; i++)
		{
			neighboring_ids.clear();
			particle_id = branches_[1]->inner_particles_.front() + i;
			if (i == 0)
			{
				neighboring_ids.push_back(branches_[0]->inner_particles_.front());
				neighboring_ids.push_back(particle_id + 1);
				neighboring_ids.push_back(particle_id + 2);
			}
			else if (i == 1)
			{
				neighboring_ids.push_back(branches_[0]->inner_particles_.front());
				neighboring_ids.push_back(particle_id - 1);
				neighboring_ids.push_back(particle_id + 1);
				neighboring_ids.push_back(particle_id + 2);
			}
			else if (2 <= i && i <= (num_ele - 3))
			{
				neighboring_ids.push_back(particle_id - 1);
				neighboring_ids.push_back(particle_id - 2);
				neighboring_ids.push_back(particle_id + 1);
				neighboring_ids.push_back(particle_id + 2);
			}
			else if (i == (num_ele - 2))
			{
				neighboring_ids.push_back(particle_id - 2);
				neighboring_ids.push_back(particle_id - 1);
				neighboring_ids.push_back(particle_id + 1);

				for (size_t k = 0; k < branches_[1]->out_edge_.size(); ++k)
				{
					size_t child_branch_id = branches_[1]->out_edge_[k];
					neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front());
				}
			}
			else if (i == (num_ele - 1))
			{
				neighboring_ids.push_back(particle_id - 1);
				neighboring_ids.push_back(particle_id - 2);

				for (size_t k = 0; k < branches_[1]->out_edge_.size(); ++k)
				{
					size_t child_branch_id = branches_[1]->out_edge_[k];
					neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front());
					neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front() + 1);
				}
			}

			for (size_t n = 0; n != neighboring_ids.size(); ++n)
			{
				Vecd displacement = pos_[particle_id] - pos_[neighboring_ids[n]];
				Neighborhood &neighborhood = particle_configuration[particle_id];
				neighbor_relation_inner(neighborhood, displacement, particle_id, neighboring_ids[n]);
			}
		}
		/** Other branches. 
		 * They are may normal branch (fully growed, has child and parent) or non-fully growed branch
		 */
		for (size_t branch_idx = 2; branch_idx != branches_.size(); ++branch_idx)
		{
			num_ele = branches_[branch_idx]->inner_particles_.size();
			size_t parent_branch_id = branches_[branch_idx]->in_edge_;
			if (!branches_[branch_idx]->is_terminated_)
			{
				/** This branch is fully growed. */
				for (size_t i = 0; i != num_ele; i++)
				{
					neighboring_ids.clear();
					particle_id = branches_[branch_idx]->inner_particles_.front() + i;
					if (i == 0)
					{
						neighboring_ids.push_back(branches_[parent_branch_id]->inner_particles_.back());
						neighboring_ids.push_back(branches_[parent_branch_id]->inner_particles_.back() - 1);

						neighboring_ids.push_back(particle_id + 1);
						neighboring_ids.push_back(particle_id + 2);
					}
					else if (i == 1)
					{
						neighboring_ids.push_back(branches_[parent_branch_id]->inner_particles_.back());
						neighboring_ids.push_back(particle_id - 1);
						neighboring_ids.push_back(particle_id + 1);
						neighboring_ids.push_back(particle_id + 2);
					}
					else if (2 <= i && i <= (num_ele - 3))
					{
						neighboring_ids.push_back(particle_id - 1);
						neighboring_ids.push_back(particle_id - 2);
						neighboring_ids.push_back(particle_id + 1);
						neighboring_ids.push_back(particle_id + 2);
					}
					else if (i == (num_ele - 2))
					{
						neighboring_ids.push_back(particle_id - 2);
						neighboring_ids.push_back(particle_id - 1);
						neighboring_ids.push_back(particle_id + 1);

						for (size_t k = 0; k < branches_[branch_idx]->out_edge_.size(); ++k)
						{
							size_t child_branch_id = branches_[branch_idx]->out_edge_[k];
							neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front());
						}
					}
					else if (i == (num_ele - 1))
					{
						neighboring_ids.push_back(particle_id - 1);
						neighboring_ids.push_back(particle_id - 2);

						for (size_t k = 0; k < branches_[branch_idx]->out_edge_.size(); ++k)
						{
							size_t child_branch_id = branches_[branch_idx]->out_edge_[k];
							neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front());
							if (branches_[child_branch_id]->inner_particles_.size() >= 2)
							{
								neighboring_ids.push_back(branches_[child_branch_id]->inner_particles_.front() + 1);
							}
						}
					}

					for (size_t n = 0; n != neighboring_ids.size(); ++n)
					{
						Vecd displacement = pos_[particle_id] - pos_[neighboring_ids[n]];
						Neighborhood &neighborhood = particle_configuration[particle_id];
						neighbor_relation_inner(neighborhood, displacement, particle_id, neighboring_ids[n]);
					}
				}
			}
			else
			{
				/** This branch is not fully growed. */
				for (size_t i = 0; i != num_ele; i++)
				{
					neighboring_ids.clear();
					particle_id = branches_[branch_idx]->inner_particles_.front() + i;
					if (i == 0)
					{
						neighboring_ids.push_back(branches_[parent_branch_id]->inner_particles_.back());
						if (branches_[parent_branch_id]->inner_particles_.size() >= 2)
							neighboring_ids.push_back(branches_[parent_branch_id]->inner_particles_.back() - 1);
					}
					else if (i == 1)
					{
						neighboring_ids.push_back(branches_[parent_branch_id]->inner_particles_.back());
						neighboring_ids.push_back(particle_id - 1);
					}
					else
					{
						neighboring_ids.push_back(particle_id - 1);
						neighboring_ids.push_back(particle_id - 2);
					}

					if (i + 1 < num_ele)
						neighboring_ids.push_back(particle_id + 1);
					if (i + 2 < num_ele)
						neighboring_ids.push_back(particle_id + 2);

					for (size_t n = 0; n != neighboring_ids.size(); ++n)
					{
						Vecd displacement = pos_[particle_id] - pos_[neighboring_ids[n]];
						Neighborhood &neighborhood = particle_configuration[particle_id];
						neighbor_relation_inner(neighborhood, displacement, particle_id, neighboring_ids[n]);
					}
				}
			}
		}
	}
	//=================================================================================================//
	size_t TreeBody::BranchLocation(size_t particle_idx)
	{
		return particle_idx < base_particles_->total_real_particles_ ? branch_locations_[particle_idx] : MaxSize_t;
	}
	//=================================================================================================//
	TreeBody::Branch::Branch(TreeBody *tree)
		: Edge<size_t, IndexVector>(tree), is_terminated_(false)
	{
		tree->branches_.push_back(this);
		tree->last_branch_id_ = id_;
	}
	//=================================================================================================//
	TreeBody::Branch::Branch(size_t parent_id, TreeBody *tree)
		: Edge<size_t, IndexVector>(parent_id, tree), is_terminated_(false)
	{
		tree->branches_[parent_id]->out_edge_.push_back(id_);
		tree->branches_.push_back(this);
		tree->last_branch_id_ = id_;
	}
	//=================================================================================================//
	TreeTerminates::TreeTerminates(SPHBody &sph_body)
		: BodyPartByParticle(sph_body, "Leaves"),
		  tree_(DynamicCast<TreeBody>(this, sph_body))
	{
		for (const auto *branch : tree_.branches_)
		{
			if (branch->is_terminated_)
			{
				size_t particle_index = branch->inner_particles_.back();
				body_part_particles_.push_back(particle_index);
			}
		}
	}
	//=================================================================================================//
}
