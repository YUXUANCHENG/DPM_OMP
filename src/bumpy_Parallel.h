#ifndef BUMPY_PARALLEL_H
#define BUMPY_PARALLEL_H

#include "DPM_Parallel.h"
#include "bumpy.h"
#include "NVE.h"

class Bumpy_Parallel : public virtual Bumpy, public virtual DPM_Parallel {
public:
	// g++ does not support multiple inheritance of constructors???
	//using Bumpy::Bumpy;
	Bumpy_Parallel(int ncells, int nt, int nprint, double l, double s) :cellPacking2D::cellPacking2D(ncells, nt, nprint, l, s) {};
	Bumpy_Parallel() = default;

	using DPM_Parallel::initialize_subsystems;
	
	virtual void hopperForces(double w0, double w, double th, double g, int closed){
		bumpy_Forces();
		if (gOn){
			for (int ci = 0; ci < NCELLS; ci++) 
				cell(ci).gravityForces(g, gDire);
		}
		hopperWallForcesDP(w0,w,th,closed);
	}

	virtual void bumpy_Forces() { 
		// reset forces
		for (int ci = 0; ci < NCELLS; ci++) {
			// reset center of mass forces
			for (int d = 0; d < NDIM; d++)
			{
				cell(ci).setCForce(d, 0.0);
				cell(ci).torque = 0.0;
			}
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

#pragma omp parallel for
		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
			subsystem[i].calculateForces_insub();
		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
			subsystem[i].calculateForces_betweensub();
		for (int ci = 0; ci < NCELLS; ci++) {
			for (int d = 0; d < NDIM; d++)
				cell(ci).setCForce(d, cell(ci).cforce(d));
		}
		addUpStress();
	}

};

#endif