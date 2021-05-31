#ifndef NVE_H
#define NVE_H

#include "cellPacking2D.h"
#include <random>

class DPMNVEsimulator{
public:
	cellPacking2D * cellpointer = nullptr;
	double init_E = 0, init_U = 0, init_K = 0;
	int dof = 2;

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
			printRoutine(count, print_frequency, t);
			NVERoutine();
			count++;
		}
	}

	virtual void printRoutine(int count, int print_frequency, double t){
		cellpointer->printRoutine(count, print_frequency, t, init_E, init_U);
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

	void subspaceManager() {
		for (int i = 0; i < cellpointer->N_systems[0] * cellpointer->N_systems[1]; i++)
			cellpointer->subsystem[i].migrate_out();
		for (int i = 0; i < cellpointer->N_systems[0] * cellpointer->N_systems[1]; i++)
			cellpointer->subsystem[i].reset_cashe();
		for (int i = 0; i < cellpointer->N_systems[0] * cellpointer->N_systems[1]; i++)
			cellpointer->subsystem[i].cashe_out(0);
		for (int i = 0; i < cellpointer->N_systems[0] * cellpointer->N_systems[1]; i++)
			cellpointer->subsystem[i].cashe_out(1);
	}

};

class DPMNVEsimulatorParallel : public DPMNVEsimulator {

public:
	DPMNVEsimulatorParallel(cellPacking2D* cell):DPMNVEsimulator(cell) {
		cell->initialize_subsystems(10, 10);
	}

	virtual void NVERoutine() {

		for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
			cellpointer->cell(ci).verletPositionUpdate(cellpointer->dt);
			cellpointer->cell(ci).updateCPos();
		}

		// reset contacts before force calculation
		cellpointer->resetContacts();
		// calculate forces
		subspaceManager();

		cellpointer->calculateForces_parallel();
		// update velocities
		verletVelocityUpdate();
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
	virtual void printRoutine(int count, int print_frequency, double t){
		if (count % print_frequency == 0)
			cout << "t = " << t << endl;
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
	virtual void printRoutine(int count, int print_frequency, double t){
		if (count % print_frequency == 0)
			cout << "t = " << t << endl;
	}
};

class DPMhopperSimulator {
public:
	cellPacking2D* cellpointer = nullptr;
	int closed = 1;
	int N_inside = 0;
	int clog_count = 0;

	DPMhopperSimulator(cellPacking2D* cell) {
		cellpointer = cell;
	}

	DPMhopperSimulator() = default;

	int hopperFlow(double w0, double w, double th, double g, double b) {
		int result;
		for (int t = 0; t < cellpointer->NT; t++) {
			if (closed == 1 && t > cellpointer->NT / 500) closed = 0;
			cellpointer->printRoutine(t, cellpointer->NPRINT, t, N_inside, closed);
			hopperRoutine(w0, w, th, g, b);
			result = checkTermination();
			if (result < 2)
				return result;
		}
		return 1;
	}

	int checkTermination() {
		if (N_inside == 0)
		{
			return 0;
		}
		if (Ke() < 1e-16 * N_inside && closed == 0)
		{
			clog_count++;
			if (clog_count > 100)
				return 1;
		}
		else
			clog_count = 0;
		return 2;
	}

	virtual double Ke() {
		return cellpointer->totalKineticEnergy();
	}

	virtual void hopperRoutine(double w0, double w, double th, double g, double b) {
		N_inside = 0;
		for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
			if (cellpointer->cell(ci).inside_hopper)
			{
				N_inside++;
				cellpointer->cell(ci).verletPositionUpdate(cellpointer->dt);
				cellpointer->cell(ci).updateCPos();
				// if still inside hopper
				if (cellpointer->cell(ci).cpos(0) > cellpointer->L.at(0) * 1.5)
					cellpointer->cell(ci).inside_hopper = 0;
			}
		}
		// reset contacts before force calculation
		cellpointer->resetContacts();
		// calculate forces
		cellpointer->hopperForces(w0, w, th, g, closed);
		// update velocities
		for (int ci = 0; ci < cellpointer->NCELLS; ci++)
			cellpointer->cell(ci).verletVelocityUpdate(cellpointer->dt, b);
	}
};
class BumpyHopperSimulator : public DPMhopperSimulator {
public:

	BumpyHopperSimulator(cellPacking2D* cell) {
		cellpointer = cell;
	}

	BumpyHopperSimulator() = default;

	virtual double Ke() {
		return cellpointer->totalKineticEnergy() + cellpointer->totalRotaionalK();
	}
	virtual void hopperRoutine(double w0, double w, double th, double g, double b) {
		N_inside = 0;
		for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
			if (cellpointer->cell(ci).inside_hopper)
				N_inside++;
		}
		cellpointer->spPosVerlet();
		cellpointer->bumpyRotation();
		for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
			if (cellpointer->cell(ci).cpos(0) > cellpointer->L.at(0) * 1.5)
				cellpointer->cell(ci).inside_hopper = 0;
		}
		// reset contacts before force calculation
		cellpointer->resetContacts();
		// calculate forces
		cellpointer->hopperForces(w0, w, th, g, closed);
		// update velocities
		cellpointer->sp_VelVerlet(b);
		cellpointer->bumpy_angularV(b);

	}
};

#endif
