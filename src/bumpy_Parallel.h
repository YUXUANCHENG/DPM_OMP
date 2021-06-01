#ifndef BUMPY_PARALLEL_H
#define BUMPY_PARALLEL_H

#include "DPM_Parallel.h"
#include "bumpy.h"
#include "NVE.h"

class Bumpy_Parallel : public Bumpy, public DPM_Parallel {
public:
	using Bumpy::Bumpy;

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

		subspaceManager();

#pragma omp for
		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
			subsystem[i].calculateForces_insub();
		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
			subsystem[i].calculateForces_betweensub();

	}

};

#endif