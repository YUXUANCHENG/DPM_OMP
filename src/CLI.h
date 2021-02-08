#ifndef CLI_H
#define CLI_H

// include files
#include "cellPacking2D.h"
#include "bumpy.h"
#include "deformableParticles2D.h"
#include <sstream>
#include <cmath>
#include <omp.h>

class DPM_CLI {
public:
	// define PI
	const double PI = 4.0 * atan(1);

	// length paramaters
	const int NT = 1e7;
	const int NPRINT = 10000;

	// simulation constants
	const double sizedev = 0.1;			        // std dev of cell sizes
	double timeStepMag = 0.005;		// time step in MD units (zeta * lenscale / forcescale)

	// disk constants
	double phiDisk = 0.5;			// initial packing fraction of disks

	// compression constants
	const double phiTarget = 1.03;			// cell packing fraction (regardless of final pressure)
	const double deltaPhi = 0.001;		// compression step size

	// force parameters
	double kl = 0.1;				// perimeter force constant
	double ka = 1.0;				// area force constant
	double gam = 0.0;				// surface tension force constant
	double kb = 0.0;				// bending energy constant
	//const double kb				= 0.1;				// bending energy constant

	const double kint = 1.0;				// interaction energy constant
	const double del = 1.0;				// width of vertices in units of l0, vertex sep on regular polygon
	const double aInitial = 0.0;				// attraction parameter to start

	// ratio of preferred perimeter^2 to preferred area
	double calA0 = 1.0;

	// tolerances
	double Ftolerance = 1e-10;			// force tolerance (for FIRE min)
	double Ptolerance = 1e-8;
	double Ktolerance = 1e-16;
	string extend;

	// system size
	int NCELLS = 16;
	int NV = 16;
	int seed = 1;
	double Lini = 1.0;

	// activity
	double T = 1000000.0;
	int frames = 50000;
	double Dr;
	double vtau = 1e-2;
	double t_scale = 1.00;

	cellPacking2D* particles;

	template <class Ptype = cellPacking2D>
	void _createParticles(char const* argv[])
	{
		int index_i;
		string index_str = argv[1];
		stringstream indexss(index_str);
		indexss >> index_i;

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");
		//double ratio = 100.0;
		kb = 0.00001 * pow(index_i + 1, 2);
		//double kl = ratio * kb;
		phiDisk = 0.78;

		// output files
		string extend = "_jammed_" + to_string(index_i) + ".txt";
		string energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF;
		produceFileName(extend, energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF);
		// instantiate object
		cout << "	** Cell packing, NCELLS = " << NCELLS << endl;
		particles = new Ptype(NCELLS, NT, NPRINT, Lini, seed);
		particles->openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

	}

	void produceFileName(string extend, string& energyF, string& jammingF, string& lengthscaleF,
		string& phiF, string& calAF, string& contactF, string& vF) {
		//string positionF = "position" + extend;
		energyF = "energy" + extend;
		jammingF = "jam" + extend;
		lengthscaleF = "length" + extend;
		phiF = "phi" + extend;
		calAF = "calA" + extend;
		contactF = "contact" + extend;
		vF = "v" + extend;
	}

	virtual void prepareSystem() {
		// Initialze the system as disks
		cout << "	** Initializing at phiDisk = " << phiDisk << endl;
		particles->initializeGel(NV, phiDisk, sizedev, del);
		particles->forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
		particles->vertexDPMTimeScale(timeStepMag);
		particles->qsIsoCompression(phiDisk, deltaPhi, Ftolerance);
	}

	void findJamming(char const* argv[]) {
		createParticles(argv);
		prepareSystem();
		particles->findJamming(deltaPhi, Ftolerance, Ptolerance);
	}

	virtual void createParticles(char const* argv[])
	{
		_createParticles<cellPacking2D>(argv);
	}
};

class Bumpy_CLI : public DPM_CLI {
public:

	virtual void createParticles(char const* argv[])
	{
		_createParticles<Bumpy>(argv);
	}

	virtual void prepareSystem() {
		// Initialze the system as disks
		cout << "	** Initializing at phiDisk = " << phiDisk << endl;
		particles->initializeGel(NV, phiDisk, sizedev, del);
		particles->vertexDPMTimeScale(timeStepMag);
		//particles->qsIsoCompression(phiDisk, deltaPhi, Ftolerance);
	}
};

#endif