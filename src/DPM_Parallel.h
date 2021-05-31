#ifndef DPM_PARALLEL_H
#define DPM_PARALLEL_H

#include "cellPacking2D.h"
#include "NVE.h"

class DPM_Parallel : public cellPacking2D {
public:
	using cellPacking2D::cellPacking2D;

	void NVEsimulation(double T, double v0, double t_scale, int frames) {
		LangevinSimulation(4000, v0, t_scale, 20);
		DPMNVEsimulatorParallel simulator = DPMNVEsimulatorParallel(this);
		simulator.NVEsimulationNoInjection(T, v0, t_scale, frames);
	}
};

#endif