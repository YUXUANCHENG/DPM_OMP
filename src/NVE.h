#ifndef NVE_H
#define NVE_H

#include "cellPacking2D.h"
#include <random>

class DPMNVEsimulator{
public:
	cellPacking2D * cellpointer = nullptr;
	double init_E = 0, init_U = 0, init_K = 0;
	int dof;

	DPMNVEsimulator(cellPacking2D* cell) {
		cellpointer = cell;
		dof = 2;
	}

	DPMNVEsimulator() = default;

	virtual void NVERoutine() {

		for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
			cellpointer->cell(ci).verletPositionUpdate(cellpointer->dt);
			cellpointer->cell(ci).updateCPos();
		}

		// reset contacts before force calculation
		cellpointer->resetContacts();
		// calculate forces
		cellpointer->calculateForces();
		// update velocities
		verletVelocityUpdate();
	}

	virtual void verletVelocityUpdate(){
		for (int ci = 0; ci < cellpointer->NCELLS; ci++)
			cellpointer->cell(ci).verletVelocityUpdate(cellpointer->dt);
	}

	virtual void NVEsimulation(double T, double v0, double t_scale, int frames) {
		injectT(v0);
		NVEsimulationNoInjection(T, v0, t_scale, frames);
	}

	virtual void NVEsimulationNoInjection(double T, double v0, double t_scale, int frames) {
		int count = 0;
		int print_frequency = floor(T / (cellpointer->dt0 * t_scale * frames));

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
		totalDof = 2 * dof * totalNV / cellpointer->NCELLS;
		double scaled_v = cellpointer->scale_v(v0);
		init_K = cellpointer->cal_temp(scaled_v);
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
	BumpyNVEsimulator() = default;

	virtual void NVERoutine() {
		// update positions
		cellpointer->spPosVerlet();
		cellpointer->bumpyRotation();
		// reset contacts before force calculation
		cellpointer->resetContacts();
		// calculate forces
		cellpointer->bumpy_Forces();
		// update velocities
		verletVelocityUpdate();
	}

	virtual void verletVelocityUpdate(){
		cellpointer->sp_VelVerlet();
		cellpointer->bumpy_angularV();
	}

	virtual void injectT(double v0) {
		double totalDof;
		totalDof = dof;
		double scaled_v = cellpointer->scale_v(v0);
		init_K = cellpointer->cal_temp(scaled_v);
		init_U = cellpointer->totalPotentialEnergy();
		init_E = totalDof * init_K + init_U;
		// Reset velocity
		cellpointer->resetV();
		cellpointer->rescal_V(init_E);
	}
};

class BumpyLangevin : public BumpyNVEsimulator {
public:
	// Mersenne twister PRNG, initialized with seed from previous random device instance
	std::mt19937 gen;
	std::normal_distribution<double> dist;

	BumpyLangevin(cellPacking2D* cell){
		cellpointer = cell;
		dof = 3;
		std::random_device rd;
		gen = std::mt19937(rd());
		dist = std::normal_distribution<double>(0, 1);
	}
	BumpyLangevin() = default;

	virtual void verletVelocityUpdate(){
		cellpointer->sp_VelVerlet_Langevin(1e-2, 2* init_K/ cellpointer->NCELLS, dist, gen);
		cellpointer->bumpy_angularV();
	}
};

class DPMLangevin : public DPMNVEsimulator {
public:
	// Mersenne twister PRNG, initialized with seed from previous random device instance
	std::mt19937 gen;
	std::normal_distribution<double> dist;

	DPMLangevin(cellPacking2D* cell){
		cellpointer = cell;
		dof = 2;
		std::random_device rd;
		gen = std::mt19937(rd());
		dist = std::normal_distribution<double>(0, 1);
	}
	DPMLangevin() = default;

	virtual void verletVelocityUpdate(){
		for (int ci = 0; ci < cellpointer->NCELLS; ci++)
			cellpointer->cell(ci).velVerlet_Langevin(cellpointer->dt, 1e-2, 2* init_K/ cellpointer->NCELLS, dist, gen);
		cellpointer->conserve_momentum();
	}
};

#endif
