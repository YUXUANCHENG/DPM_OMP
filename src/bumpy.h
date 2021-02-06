#ifndef BUMPY_H
#define BUMPY_H

#include "cellPacking2D.h"
#include "NVE.h"

class Bumpy : public cellPacking2D {
public:
	using cellPacking2D::cellPacking2D;
	virtual void NVEsimulation(double T, double v0, double t_scale, int frames);
	void initializeGel(int NV, double phiDisk, double sizeDispersion, double delval);
	virtual void printRoutine(int count, int print_frequency, double t, double init_E, double init_U);
	virtual void resetV();
	virtual void fireMinimizeF();
};

#endif