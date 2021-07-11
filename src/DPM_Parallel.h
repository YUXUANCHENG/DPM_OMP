#ifndef DPM_PARALLEL_H
#define DPM_PARALLEL_H

#include "cellPacking2D.h"
#include "NVE.h"

class DPM_Parallel : public virtual cellPacking2D {
public:
	using cellPacking2D::cellPacking2D;

	virtual void initialize_subsystems(int N_x, int N_y);
	void reset_subsystems();
	void delete_subsystems();
	void split_into_subspace();
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
			cell(ci).shapeForces();
		}
		addUpStress();
	}

	virtual void hopperForces(double w0, double w, double th, double g, int closed){
		calculateForces();
		if (gOn){
#pragma omp parallel for
			for (int ci = 0; ci < NCELLS; ci++) 
				cell(ci).gravityForces(g, gDire);
		}
		hopperWallForcesDP(w0,w,th,closed);
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

};

#endif