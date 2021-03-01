#ifndef BUMPY_H
#define BUMPY_H

#include "cellPacking2D.h"
#include "NVE.h"

class Bumpy : public cellPacking2D {
public:
	using cellPacking2D::cellPacking2D;
	void printV();
	virtual void compressToInitial(double phiTarget, double deltaPhi, double Ftol);
	virtual void qsIsoCompression(double phiDisk, double deltaPhi, double Ftolerance);
	virtual void NVEsimulation(double T, double v0, double t_scale, int frames);
	virtual void LangevinSimulation(double T, double v0, double t_scale, int frames);
	virtual void printRoutine(int count, int print_frequency, double t, double init_E, double init_U);
	virtual void resetV();
	virtual void fireMinimizeF(double Ftol, double& Ftest, double& Ktest);
	virtual int hopperSimulation(double w0, double w, double th, double g, double b);
	virtual void hopperForces(double w0, double w, double th, double g, int closed);
};

#endif