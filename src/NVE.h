#ifndef NVE_H
#define NVE_H

#include "cellPacking2D.h"

class DPMNVEsimulator{
public:
	cellPacking2D * cellpointer = nullptr;
	double init_E = 0, init_U = 0;
	int dof = 2;

	DPMNVEsimulator(cellPacking2D* cell) {
		cellpointer = cell;
		dof = 4;
	}

	DPMNVEsimulator() = default;

	virtual void NVERoutine() {
		int ci;

		for (ci = 0; ci < cellpointer->NCELLS; ci++) {
			cellpointer->cell(ci).verletPositionUpdate(cellpointer->dt);
			cellpointer->cell(ci).updateCPos();
		}

		// reset contacts before force calculation
		cellpointer->resetContacts();

		// calculate forces
		cellpointer->calculateForces();

		// update velocities
		for (ci = 0; ci < cellpointer->NCELLS; ci++)
			cellpointer->cell(ci).verletVelocityUpdate(cellpointer->dt);
	}

	virtual void NVEsimulation(double T, double v0, double t_scale, int frames) {
		int ci, vi, d;
		int count = 0;
		double U, K, rv;
		int print_frequency = floor(T / (cellpointer->dt0 * t_scale * frames));

		injectT(v0);

		for (double t = 0.0; t < T; t = t + cellpointer->dt0 * t_scale) {
			cellpointer->printRoutine(count, print_frequency, t, init_E, init_U);
			NVERoutine();
			count++;
		}
	}

	virtual void injectT(double v0) {
		int ci;
		int totalNV = 0;
		double totalDof;
		for (ci = 0; ci < cellpointer->NCELLS; ci++) {
			//system_mass += cell(ci).getNV() * PI * pow(0.5 * cell(ci).getdel() * cell(ci).getl0(), 2);
			totalNV += cellpointer->cell(ci).getNV();
		}
		totalDof = dof * totalNV / cellpointer->NCELLS;
		double scaled_v = cellpointer->scale_v(v0);
		double init_K = cellpointer->cal_temp(scaled_v);
		init_U = cellpointer->totalPotentialEnergy();
		init_E = totalDof * init_K + init_U;
		// Reset velocity
		cellpointer->resetV();
		cellpointer->rescal_V(init_E);
	}

};

class BumpyNVEsimulator : public DPMNVEsimulator {
public:

	BumpyNVEsimulator(cellPacking2D* cell){
		cellpointer = cell;
		dof = 3;
	}

	virtual void NVERoutine() {
		// update positions
		cellpointer->spPosVerlet();
		cellpointer->bumpyRotation();

		// reset contacts before force calculation
		cellpointer->resetContacts();

		// calculate forces
		cellpointer->bumpy_Forces();
		//calculateForces();

		// update velocities
		cellpointer->sp_VelVerlet();
		cellpointer->bumpy_angularV();
	}

	virtual void injectT(double v0) {
		double totalDof;
		totalDof = dof;
		double scaled_v = cellpointer->scale_v(v0);
		double init_K = cellpointer->cal_temp(scaled_v);
		init_U = cellpointer->totalPotentialEnergy();
		init_E = totalDof * init_K + init_U;
		// Reset velocity
		cellpointer->resetV();
		cellpointer->rescal_V(init_E);
	}
};

#endif
