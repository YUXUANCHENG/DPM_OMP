#ifndef DPM_PARALLEL_H
#define DPM_PARALLEL_H

#include "cellPacking2D.h"
#include "NVE.h"
extern bool constPressureFlag;
extern bool pinFlag;
class DPM_Parallel : public virtual cellPacking2D {
public:
	using cellPacking2D::cellPacking2D;

	virtual void initialize_subsystems(int N_x, int N_y);
	void reset_subsystems();
	void delete_subsystems();
	virtual void split_into_subspace();
	void cashe_into(int i, vector<cvpair*>& cash_list);
	void migrate_into(int i, cvpair* const& migration);
	int look_for_new_box(cvpair* pair);
	double transformPos(cvpair* pair, int direction);

/*
	void NVEsimulation(double T, double v0, double t_scale, int frames) {
		LangevinSimulation(4000, v0, t_scale, 20);
		DPMNVEsimulatorParallel simulator = DPMNVEsimulatorParallel(this);
		simulator.NVEsimulationNoInjection(T, v0, t_scale, frames);
	}
*/

	virtual void scaleLengths(double val) {
		cellPacking2D::scaleLengths(val);
		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
			subsystem[i].cal_cashed_fraction();
	}

	virtual void initialize_subsystems() {
		initialize_subsystems(N_systems.at(0) , N_systems.at(1));
	}


	void subspaceManager() {
		for (int i = 0; i < N_systems[0] * N_systems[1]; i++)
			subsystem[i].migrate_out();
		for (int i = 0; i < N_systems[0] * N_systems[1]; i++)
			subsystem[i].reset_cashe();
		for (int i = 0; i < N_systems[0] * N_systems[1]; i++)
			subsystem[i].cashe_out(0);
		for (int i = 0; i < N_systems[0] * N_systems[1]; i++)
			subsystem[i].cashe_out(1);
	}

	virtual void calculateForces() {

		// reset forces
		for (int ci = 0; ci < NCELLS; ci++) {
			// reset center of mass forces
			for (int d = 0; d < NDIM; d++)
				cell(ci).setCForce(d, 0.0);

			// reset vertex forces and interaction energy
			for (int vi = 0; vi < cell(ci).getNV(); vi++) {
				// forces
				for (int d = 0; d < NDIM; d++)
					cell(ci).setVForce(vi, d, 0.0);

				// energies
				cell(ci).setUInt(vi, 0.0);
			}
		}

		resetContacts();

		subspaceManager();
#pragma omp parallel for schedule (dynamic, 4)
		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
			subsystem[i].calculateForces_insub();
#pragma omp parallel for schedule (dynamic, 4)
		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
			subsystem[i].calculateForces_betweensub();
#pragma omp parallel for
		for (int ci = 0; ci < NCELLS; ci++) {
			if (cell(ci).inside_hopper)
				cell(ci).shapeForces();
		}
		addUpStress();
	}

	int countTotalVertices()
	{
		int count = 0;
		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
			count += subsystem[i].resident_cells.size();
		return count;
	}

	virtual void hopperForces(double w0, double w, double th, double g, int closed){
		calculateForces();
		if (gOn){
#pragma omp parallel for
			for (int ci = 0; ci < NCELLS; ci++) 
				if (cell(ci).inside_hopper)
					cell(ci).gravityForces(g, gDire);
		}
		if (!pinFlag)
			hopperWallForcesDP(w0,w,th,closed);
		else
			hopperWallForcesPining(w0,w,th,closed);
		if (constPressureFlag)
			pistonForce(w0, w, th, g);
	}


	void addUpStress()
	{
		// reset virial stresses to 0
		sigmaXX = 0.0;
		sigmaXY = 0.0;
		sigmaYX = 0.0;
		sigmaYY = 0.0;

		// reset contacts before force calculation
		//resetContacts();
		Ncc = 0;
		Nvv = 0;

		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
		{
			sigmaXX += subsystem[i].sigmaXX;
			sigmaXY += subsystem[i].sigmaXY;
			sigmaYX += subsystem[i].sigmaYX;
			sigmaYY += subsystem[i].sigmaYY;

			Nvv += subsystem[i].Nvv;
			// need to fix later
			Ncc += subsystem[i].Ncc;
		}
	}

	virtual void printRoutine(int count, int print_frequency, double t, double init_E, double init_U)
	{
		cellPacking2D::printRoutine(count, print_frequency, t, init_E, init_U);
		if (count % print_frequency == 0)
			cout << "total vetex count = " << countTotalVertices() << endl;
	}

};


class DPM_Parallel_frictionless : public virtual DPM_Parallel {
public:
	using DPM_Parallel::DPM_Parallel;
	virtual void split_into_subspace();

	virtual void calculateForces() {
				// reset forces
		for (int ci = 0; ci < NCELLS; ci++) {
			// reset center of mass forces
			for (int d = 0; d < NDIM; d++)
				cell(ci).setCForce(d, 0.0);

			// reset vertex forces and interaction energy
			for (int vi = 0; vi < cell(ci).getNV(); vi++) {
				// forces
				for (int d = 0; d < NDIM; d++)
					cell(ci).setVForce(vi, d, 0.0);

				// energies
				cell(ci).setUInt(vi, 0.0);
				cell(ci).vertexEdgeContact[vi] = -1;
			}
		}

		resetContacts();

		subspaceManager();
// #pragma omp parallel for schedule (dynamic, 4)
// 		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
// 			subsystem[i].calculateEdgeForces_insub();
// #pragma omp parallel for schedule (dynamic, 4)
// 		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
// 			subsystem[i].calculateEdgeForces_betweensub();

#pragma omp parallel for schedule (dynamic, 4)
		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
			subsystem[i].calculateForces_insub();
#pragma omp parallel for schedule (dynamic, 4)
		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
			subsystem[i].calculateForces_betweensub();
#pragma omp parallel for
		for (int ci = 0; ci < NCELLS; ci++) {
			if (cell(ci).inside_hopper)
				cell(ci).shapeForces();
		}
		addUpStress();
		// resolveForces();
	}
	
	void resolveForces()
	{
#pragma omp parallel for schedule (dynamic, 4)
		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
			subsystem[i].calculateEdgeForces_insub();
#pragma omp parallel for schedule (dynamic, 4)
		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
			subsystem[i].calculateEdgeForces_betweensub();
// 		for (auto const& x : collisionMap)
// 		{
// 			const cvpair & onTheLeft = x.first;
// 			const std::vector<cvpair> & onTheRightList = x.second;
// 			cvpair smallestDistPair;
// 			double smallestDist = 1e9;
// 			for (cvpair onTheRight : onTheRightList)
// 			{
// 				// double vertexDist = abs(subsystem[onTheLeft.boxid].vertexEdgeDist(&onTheLeft, &onTheRight));
// 				// if (vertexDist < smallestDist)
// 				{
// 					smallestDistPair = onTheRight;
// 					// smallestDist = vertexDist;
// 				}
// 			// }
// 			vector<VECTOR2> v, e;
// 			v.resize(3);
// 			VECTOR2 v0, v1, vColission; 
// 			int nextVi = (smallestDistPair.vi + 1) % cell(smallestDistPair.ci).getNV();
// 			for (int d = 0; d < 2; d++)
// 			{
// 				vColission[d] = cell(onTheLeft.ci).vpos(onTheLeft.vi,d);
// 				v0[d] = cell(smallestDistPair.ci).vpos(smallestDistPair.vi,d);
// 				v1[d] = cell(smallestDistPair.ci).vpos(nextVi,d);
// 			}
// 			// v[0][0] = pointer_to_system->cell(onTheLeft->ci).vpos(onTheLeft->vi,0);
// 			// v[0][1] = pointer_to_system->cell(onTheLeft->ci).vpos(onTheLeft->vi,1);
// 			// v[1][0] = pointer_to_system->cell(onTheRight->ci).vpos(onTheRight->vi,0);
// 			// v[1][1] = pointer_to_system->cell(onTheRight->ci).vpos(onTheRight->vi,1);
// 			// int nextVi = (onTheRight->vi + 1) % pointer_to_system->cell(onTheRight->ci).getNV();
// 			// v[2][0] = pointer_to_system->cell(onTheRight->ci).vpos(nextVi,0);
// 			// v[2][1] = pointer_to_system->cell(onTheRight->ci).vpos(nextVi,1);

// 			v[0] = vColission;
// 			v[1] = v0;
// 			v[2] = v1;

// 			e.resize(2);
// 			e[0] = v[2] - v[1];
// 			e[1] = v[0] - v[1];
// 			deformableParticles2D& leftCell = cell(onTheLeft.ci);
// 			deformableParticles2D& rightCell = cell(smallestDistPair.ci);
// 			double eps = 0.5 * (leftCell.del * leftCell.l0 + rightCell.del * rightCell.l0) * cutoff;
// 			VECTOR6 forces = -1 * gradient(v, e, eps);
// #pragma omp critical
// {
// 			for (int d = 0; d < 2; d++)
// 			{
// 				leftCell.setVForce(onTheLeft.vi, d, leftCell.vforce(onTheLeft.vi, d) + forces[d]);
// 				rightCell.setVForce(smallestDistPair.vi, d, rightCell.vforce(smallestDistPair.vi, d) + forces[d + 2]);
// 				rightCell.setVForce(nextVi, d, rightCell.vforce(nextVi, d) + forces[d + 4]);
// 			}
// }
// 			}
// 			collisionMap[onTheLeft].clear();
// 			collisionMap[onTheLeft].resize(0);
// 		}
// 		collisionMap.clear();
	}

	VECTOR6 springLengthGradient(const vector<VECTOR2>& v,
														const vector<VECTOR2>& e,
														const VECTOR2& n)
	{
	const MATRIX2x6 nPartial = normalGradientVF(e);
	const VECTOR2 tvf = v[0] - v[2];

	MATRIX2x6 tvfPartial;
	tvfPartial.setZero();
	tvfPartial(0,0) = tvfPartial(1,1) = 1.0;
	tvfPartial(0,4) = tvfPartial(1,5) = -1.0;

	//f = nPartial' * (v2 - v0) + tvfPartial' * n;
	return nPartial.transpose() * tvf + tvfPartial.transpose() * n;
	}

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	virtual VECTOR6 gradient(vector<VECTOR2> & v, vector<VECTOR2> & e, double eps, int mode)
	{
	// convert to vertices and edges
	//   vector<VECTOR2> v;
	//   vector<VECTOR2> e;
	//   getVerticesAndEdges(x, v, e);
	
	// get the normal
	VECTOR2 n = VECTOR2(e[0][1],-e[0][0]);
	n = n / n.norm();
	
	// get the spring length, non-zero rest-length
	const VECTOR2 tvf = v[0] - v[1];
	double normalDist = n.dot(tvf);
	double springLength = normalDist - eps;
	if (mode && normalDist < 0)
		springLength = normalDist + eps;
	return 2.0 * springLength * springLengthGradient(v,e,n);
	}

	///////////////////////////////////////////////////////////////////////
	// gradient of the triangle normal, vertex-face case
	///////////////////////////////////////////////////////////////////////
	MATRIX2x6 normalGradientVF(const std::vector<VECTOR2>& e)
	{
	//crossed = cross(e2, e0);
	VECTOR2 crossed = VECTOR2(e[0][1],-e[0][0]);
	double crossNorm = crossed.norm();
	const double crossNormCubedInv = 1.0 / pow(crossed.dot(crossed), 1.5);
	MATRIX2x6 crossMatrix = crossGradientVF(e);

	//final = zeros(3,12);
	//for i = 1:12
	//  crossColumn = crossMatrix(:,i);
	//  final(:,i) = (1 / crossNorm) * crossColumn - ((crossed' * crossColumn) / crossNormCubed) * crossed;
	//end
	MATRIX2x6 result;
	for (int i = 0; i < 6; i++)
	{
		const VECTOR2 crossColumn = crossMatrix.col(i);
		result.col(i) = (1.0 / crossNorm) * crossColumn - 
						((crossed.dot(crossColumn)) * crossNormCubedInv) * crossed;
	}
	return result;
	}

	MATRIX2x6 crossGradientVF(const std::vector<VECTOR2>& e)
	{
		MATRIX2x6 crossMatrix;

		crossMatrix.col(0) = VECTOR2(0, 0); 
		crossMatrix.col(1) = VECTOR2(0, 0); 
		crossMatrix.col(2) = VECTOR2(0, 1); 
		crossMatrix.col(3) = VECTOR2(-1, 0);
		crossMatrix.col(4) = VECTOR2(0, -1);
		crossMatrix.col(5) = VECTOR2(1, 0);
		return crossMatrix;
	}
};

#endif