#ifndef jamming_h
#define jamming_h


// include files
#include "cellPacking2D.h"
#include "deformableParticles2D.h"
#include <sstream>


// use std name space
using namespace std;

class jamming {
private:

	// define PI
	const double PI = 4.0 * atan(1);

	// length paramaters
	const int NT = 1e7;
	const int NPRINT = 10000;

	// simulation constants
	const double sizedev = 0.1;			        // std dev of cell sizes
	const double timeStepMag = 0.005;		// time step in MD units (zeta * lenscale / forcescale)

	// disk constants
	//const double phiDisk	 		= 0.65;			// initial packing fraction of disks
	const double phiDisk = 0.6;

	// compression constants
	const double phiTarget = 1.03;			// cell packing fraction (regardless of final pressure)
	const double deltaPhi = 0.001;		// compression step size

	// gelation constants
	const double phiGel = 0.3;			// final packing fraction
	const double gelRate = 1e-4;			// rate of size decrease (i.e. area loss relative to initial box area)
	const double varPerimRate = 0.01;			// rate of relaxation to deformed perimeter
	const double aGelation = 0.05;			// attraction parameter during gelation sim

	// force parameters
	double kl = 1.0;				// perimeter force constant
	double ka = 1.0;				// area force constant
	double gam = 0.0;				// surface tension force constant
	double kb = 0.0;				// bending energy constant
	//const double kb				= 0.1;				// bending energy constant

	const double kint = 1.0;				// interaction energy constant
	const double del = 1.0;				// width of vertices in units of l0, vertex sep on regular polygon
	const double aInitial = 0.0;				// attraction parameter to start


	// ratio of preferred perimeter^2 to preferred area
	const double calA0 = 1.1;

	// tolerances
	const double Ftolerance = 1e-8;			// force tolerance (for FIRE min)
	const double Ptolerance = 1e-8;
	const double Ktolerance = 1e-16;
	string extend;


public:
	
	void unjam()
	{

		// output files
		string positionF = "position.txt";
		string energyF = "energy.txt";
		string jammingF = "jam.txt";
		string lengthscaleF = "length.txt";
		string phiF = "phi.txt";
		string calAF = "calA.txt";
		string contactF = "contact.txt";
		string vF = "v.txt";

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");

		// system size
		int NCELLS = 16;
		int NV = 16;
		int seed = 5;
		double Lini = 1.0;

		// activity
		double T = 100.0;
		double v0 = 1.0;
		double Dr = 10.0;
		double vtau = 1e-2;
		double t_scale = 1.00;

		// instantiate object
		cout << "	** Cell packing, NCELLS = " << NCELLS << endl;
		cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

		// open position output file
		cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

		// Initialze the system as disks
		cout << "	** Initializing at phiDisk = " << phiDisk << endl;
		cell_group.initializeGel(NV, phiDisk, sizedev, del);

		// change to DPM
		cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

		// set time scale
		cell_group.vertexDPMTimeScale(timeStepMag);



		// Compress then relax by FIRE
		cout << " Compress then relax by FIRE " << endl;

		double phiTargetTmp = 0.8;
		double deltaPhiTmp = 0.001;

		//cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp);
		cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);

		cellPacking2D jammed_state;
		cell_group.saveState(jammed_state);

		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 10; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;
				v0 = 0.1 + double(i) * 0.1;;
				Dr = 1.0 + double(j) * 1.0;
				v0PrintObject << v0 << "," << Dr << endl;

				extend = "_" + to_string(i) + to_string(j) + ".txt";

				// output files
				//string positionF = "position" + extend;
				string energyF = "energy" + extend;
				string jammingF = "jam" + extend;
				string lengthscaleF = "length" + extend;
				string phiF = "phi" + extend;
				string calAF = "calA" + extend;
				string contactF = "contact" + extend;
				string vF = "v" + extend;

				cell_group.loadState(jammed_state);
				cell_group.closeF();
				cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

				cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
				cell_group.activityCOM(T, v0, Dr, vtau, t_scale);
				cell_group.relaxF(Ktolerance, Ftolerance, Ptolerance);
			}
		}

		cout << "	** FINISHED **   " << endl;
	};


	void unjam_decompressed(double d_phi)
	{

		// output files
		string positionF = "position.txt";
		string energyF = "energy.txt";
		string jammingF = "jam.txt";
		string lengthscaleF = "length.txt";
		string phiF = "phi.txt";
		string calAF = "calA.txt";
		string contactF = "contact.txt";
		string vF = "v.txt";

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");

		// system size
		int NCELLS = 16;
		int NV = 16;
		int seed = 5;
		double Lini = 1.0;

		// activity
		double T = 100.0;
		double v0 = 1.0;
		double Dr = 10.0;
		double vtau = 1e-2;
		double t_scale = 1.00;

		// instantiate object
		cout << "	** Cell packing, NCELLS = " << NCELLS << endl;
		cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

		// open position output file
		cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

		// Initialze the system as disks
		cout << "	** Initializing at phiDisk = " << phiDisk << endl;
		cell_group.initializeGel(NV, phiDisk, sizedev, del);

		// change to DPM
		cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

		// set time scale
		cell_group.vertexDPMTimeScale(timeStepMag);



		// Compress then relax by FIRE
		cout << " Compress then relax by FIRE " << endl;

		
		cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);
		//double phiTargetTmp = 0.8;
		double phiTargetTmp = cell_group.packingFraction() - d_phi;
		double deltaPhiTmp = 0.001;
		cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance, Ktolerance);

		cellPacking2D jammed_state;
		cell_group.saveState(jammed_state);

		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 10; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;
				v0 = 0.1 + double(i) * 0.1;
				Dr = 1.0 + double(j) * 1.0;
				v0PrintObject << v0 << "," << Dr << endl;

				extend = "_" + to_string(i) + to_string(j) + ".txt";

				// output files
				//string positionF = "position" + extend;
				string energyF = "energy" + extend;
				string jammingF = "jam" + extend;
				string lengthscaleF = "length" + extend;
				string phiF = "phi" + extend;
				string calAF = "calA" + extend;
				string contactF = "contact" + extend;
				string vF = "v" + extend;

				cell_group.loadState(jammed_state);
				cell_group.closeF();
				cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

				cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
				cell_group.activityCOM(T, v0, Dr, vtau, t_scale);
				cell_group.relaxF(Ktolerance, Ftolerance, Ptolerance);
			}
		}

		cout << "	** FINISHED **   " << endl;
	};

	void unjam_kb()
	{

		// output files
		string positionF = "position.txt";
		string energyF = "energy.txt";
		string jammingF = "jam.txt";
		string lengthscaleF = "length.txt";
		string phiF = "phi.txt";
		string calAF = "calA.txt";
		string contactF = "contact.txt";
		string vF = "v.txt";

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");

		// system size
		int NCELLS = 16;
		int NV = 16;
		int seed = 5;
		double Lini = 1.0;

		// activity
		double T = 100.0;
		double v0 = 1.0;
		double Dr = 10.0;
		double vtau = 1e-2;
		double t_scale = 1.00;

		for (int i = 0; i < 1; i++) {
			for (int j = 0; j < 10; j++) {


				cout << "Loop i, j = " << i << "," << j << endl;
				//v0 = 0.1 + double(i) * 0.1;;
				//Dr = 1.0 + double(j) * 1.0;
				v0 = 0.5;
				Dr = 5.0;
				kb = 0.0 + double(j) * 0.1;
				v0PrintObject << v0 << "," << Dr << "," << kb << endl;

				extend = "_" + to_string(i) + to_string(j) + ".txt";


				// instantiate object
				cout << "	** Cell packing, NCELLS = " << NCELLS << endl;
				cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

				// open position output file
				cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

				// Initialze the system as disks
				cout << "	** Initializing at phiDisk = " << phiDisk << endl;
				cell_group.initializeGel(NV, phiDisk, sizedev, del);

				// change to DPM
				cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

				// set time scale
				cell_group.vertexDPMTimeScale(timeStepMag);



				// Compress then relax by FIRE
				cout << " Compress then relax by FIRE " << endl;

				//double phiTargetTmp = 0.8;
				//double deltaPhiTmp = 0.001;

				//cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp);
				cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);

				cellPacking2D jammed_state;
				cell_group.saveState(jammed_state);

				// output files
				//string positionF = "position" + extend;
				string energyF = "energy" + extend;
				string jammingF = "jam" + extend;
				string lengthscaleF = "length" + extend;
				string phiF = "phi" + extend;
				string calAF = "calA" + extend;
				string contactF = "contact" + extend;
				string vF = "v" + extend;

				cell_group.loadState(jammed_state);
				cell_group.closeF();
				cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

				cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
				cell_group.activityCOM(T, v0, Dr, vtau, t_scale);
				cell_group.relaxF(Ktolerance, Ftolerance, Ptolerance);
			}
		}

		cout << "	** FINISHED **   " << endl;
	};
















































};


#endif