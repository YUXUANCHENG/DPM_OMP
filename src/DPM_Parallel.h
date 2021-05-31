#ifndef DPM_PARALLEL_H
#define DPM_PARALLEL_H

#include "cellPacking2D.h"
#include "NVE.h"

class DPM_Parallel : public cellPacking2D {
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
	virtual void calculateForces();
	void subspaceManager();
/*
	void NVEsimulation(double T, double v0, double t_scale, int frames) {
		LangevinSimulation(4000, v0, t_scale, 20);
		DPMNVEsimulatorParallel simulator = DPMNVEsimulatorParallel(this);
		simulator.NVEsimulationNoInjection(T, v0, t_scale, frames);
	}
*/
	void initializeGel(int NV, double phiDisk, double sizeDispersion, double delval) {
		cellPacking2D::initializeGel(NV, phiDisk, sizeDispersion, delval);
		initialize_subsystems(10, 10);
	}
};

#endif