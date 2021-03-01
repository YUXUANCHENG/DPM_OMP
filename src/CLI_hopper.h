#ifndef CLI_HOPPER_H
#define CLI_HOPPER_H

#include "CLI.h"

class DPM_Hopper_CLI : public DPM_CLI {
public:

	double smallRadius = 0.5;			// radius fo smaller particles (diameter is length unit)
	double sizeRatio = 1.4;			// ratio of small diameter to large diameter
	//const double sizeRatio = 1;
	double w0 = 20.0;			// width of hopper reservoir (in units of small diameter)

	double th = PI / 4.0;		// hopper angle (pi - th = deflection angle from horizontal)
	double phi0 = 0.4;			// initial packing fraction
	double b = 0.1;
	double g = 0.05;
	vector<double> radii;
	double w_scale = 2;
	double w;

	DPM_Hopper_CLI() {
		NT = 5e6;			// number of time steps for flow simulation
		NPRINT = 1e3;			// number of steps between printing
		kl = 10.0;
		ka = 10.0;
		kb = 10.0;
		g = 0.05;
		Lini = 0.1 * w0;
		NCELLS = 64;
		NV = 16;
		radii = vector<double>(NCELLS, 0.0);
		for (int ci = 0; ci < NCELLS; ci++) {
			if (ci % 2 == 0)
				radii.at(ci) = smallRadius;
			else
				radii.at(ci) = smallRadius * sizeRatio;
		}
	}

	virtual void setIndex(char const* argv[]) {
		string index_i_str = argv[1];
		string index_j_str = argv[2];
		stringstream index_i_ss(index_i_str);
		stringstream index_j_ss(index_j_str);
		index_i_ss >> index_i;
		index_j_ss >> index_j;
	}

	virtual void setKB() {
		;
	}

	virtual void prepareSystem() {
		w_scale = 1.5 + 0.05 * index_j;
		w = w_scale * (1 + sizeRatio) / 2;
		// Initialze the system as disks
		cout << "	** Initializing hopper " << endl;
		particles->initializeHopperDP(radii, w0, w, th, Lini, NV);
		particles->forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
		particles->vertexDPMTimeScale(timeStepMag);
		particles->closeF();
	}

	virtual void hopperFLow(char const* argv[])
	{
		qscompress(argv);
		_hopperFlow();
	}

	virtual void _hopperFlow() {
		extend = "_" + to_string(index_i) + to_string(index_j) + ".txt";
		string energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF;
		produceFileName(extend, energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF);
		particles->openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

		int result = particles->hopperSimulation(w0, w, th, g, b);
		v0PrintObject << kl << "," << gam << "," << g << "," << w_scale << "," << result << endl;
		cout << "	** FINISHED **   " << endl;
		//clogPrintObject.close();
		v0PrintObject.close();
	}
};

class Bumpy_Hopper_CLI : public DPM_Hopper_CLI {
public:
	Bumpy_Hopper_CLI() : DPM_Hopper_CLI() {
	kl = 0;
	ka = 0;
	kb = 0;
	}

	virtual void createParticles(char const* argv[])
	{
		_createParticles<Bumpy>(argv);
	}
};

#endif