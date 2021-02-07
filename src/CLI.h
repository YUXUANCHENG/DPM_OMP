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
	int frames = 1000;

	// simulation constants
	const double sizedev = 0.1;			        // std dev of cell sizes
	double timeStepMag = 0.005;		// time step in MD units (zeta * lenscale / forcescale)

	// disk constants
	//const double phiDisk	 		= 0.65;			// initial packing fraction of disks
	double phiDisk = 0.5;

	// compression constants
	const double phiTarget = 1.03;			// cell packing fraction (regardless of final pressure)
	const double deltaPhi = 0.001;		// compression step size

	// gelation constants
	const double phiGel = 0.3;			// final packing fraction
	const double gelRate = 1e-4;			// rate of size decrease (i.e. area loss relative to initial box area)
	const double varPerimRate = 0.01;			// rate of relaxation to deformed perimeter
	const double aGelation = 0.05;			// attraction parameter during gelation sim

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
	double calA0 = 1.03;

	// tolerances
	double Ftolerance = 1e-10;			// force tolerance (for FIRE min)
	double Ptolerance = 1e-8;
	double Ktolerance = 1e-16;
	string extend;

	template <class Ptype = cellPacking2D>
	void _find_jamming(char const* argv[])
	{
		int index;
		string index_str = argv[1];
		stringstream indexss(index_str);
		indexss >> index;

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");

		// system size
		int NCELLS = 16;
		int NV = 16;
		//int NV = 32;
		//int seed = 5;
		int seed = 1;
		double Lini = 1.0;

		// activity
		double T = 1000000.0;
		int frames = 50000;
		double Dr;
		double vtau = 1e-2;
		double t_scale = 1.00;

		//Dr = 1.0 + double(j) * 1.0;
		Dr = 1e-2;

		double ratio = 100.0;
		//ka = 10;

		for (int i = index; i < index + 1; i++) {

			//double calA0 = 1.12;
			double calA0 = 1.0;
			double kb = 0.00001 * pow(i + 1, 2);
			//double kl = ratio * kb;
			double kl = 0.1;

			// output files
			string extend = "_jammed_" + to_string(i) + ".txt";
			//string positionF = "position" + extend;
			string energyF = "energy" + extend;
			string jammingF = "jam" + extend;
			string lengthscaleF = "length" + extend;
			string phiF = "phi" + extend;
			string calAF = "calA" + extend;
			string contactF = "contact" + extend;
			string vF = "v" + extend;

			// instantiate object
			cout << "	** Cell packing, NCELLS = " << NCELLS << endl;
			Ptype particles(NCELLS, NT, NPRINT, Lini, seed);

			// open position output file
			particles.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
			phiDisk = 0.85;
			// Initialze the system as disks
			cout << "	** Initializing at phiDisk = " << phiDisk << endl;
			particles.initializeGel(NV, phiDisk, sizedev, del);
			//particles.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
			particles.vertexDPMTimeScale(timeStepMag);
			particles.qsIsoCompression(phiDisk, deltaPhi, Ftolerance);
			particles.findJamming(deltaPhi, Ftolerance, Ptolerance);

			Ptype jammed_state;
			particles.saveState(jammed_state);
		}
	}

	void find_jamming(char const* argv[])
	{
		_find_jamming<cellPacking2D>(argv);
	}
};

class Bumpy_CLI: public DPM_CLI {
public:
	void find_jamming(char const* argv[])
	{
		_find_jamming<Bumpy>(argv);
	}
};

#endif