#ifndef jamming_h
#define jamming_h


// include files
#include "cellPacking2D.h"
#include "deformableParticles2D.h"
#include <sstream>
#include <cmath>
#include <omp.h>

// use std name space
using namespace std;

class jamming {
private:

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
		double T = 500.0;
		double v0;
		double Dr;
		double vtau = 1e-2;
		double t_scale = 1.00;
		calA0 = 1.17;

		//timeStepMag = 0.001;

		// instantiate object
		cout << "	** Cell packing, NCELLS = " << NCELLS << endl;
		cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

		// open position output file
		cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
		phiDisk = 0.7;
		// Initialze the system as disks
		cout << "	** Initializing at phiDisk = " << phiDisk << endl;
		cell_group.initializeGel(NV, phiDisk, sizedev, del);

		// change to DPM
		cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

		// set time scale
		cell_group.vertexDPMTimeScale(timeStepMag);



		// Compress then relax by FIRE
		cout << " Compress then relax by FIRE " << endl;


		cell_group.findJamming(deltaPhi, Ftolerance, Ptolerance);

		cellPacking2D jammed_state;
		cell_group.saveState(jammed_state);

		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 10; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;

				calA0 = 1.17 + double(i) * 0.01;
				v0 = 0.04;
				//v0 = 0.01 + double(i) * 0.01;
				//Dr = 1.0 + double(j) * 1.0;
				Dr = 1e-3;
				//kb = 0.0 + double(j) * 0.005;
				kb = 0.0;
				//kl = 0.1;
				kl = 0.05 + double(j) * 0.05;
				v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << NCELLS << endl;

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
				//cell_group.relaxF(Ktolerance, Ftolerance, Ptolerance);
				cell_group.activityCOM_brownian(T, v0, Dr, vtau, t_scale, frames);
				//cell_group.relaxF(Ktolerance, Ftolerance, Ptolerance);
			}
		}

		cout << "	** FINISHED **   " << endl;
	};



	void unjam_N()
	{
		int frames = 20000;

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");

		for (int i = 5; i < 6; i++) {
			for (int j = 0; j < 2; j++) {

				// system size
				int NCELLS = pow(2, i + 5);
				int NV = 16;
				int seed = 5;
				double Lini = 1.0;

				// activity
				double T = 10000.0;
				double v0;
				double Dr;
				double vtau = 1e-2;
				double t_scale = 1.00;

				cout << "Loop i, j = " << i << "," << j << endl;

				v0 = 0.1;
				//v0 = 0.1 + double(i) * 0.1;
				//Dr = 1.0 + double(j) * 1.0;
				Dr = 1e-3;
				kl = 0.1;
				kb = 0.0 + double(j) * 0.03;
				//kb = 0.0;

				// output files
				extend = "_jammed_" + to_string(i) + to_string(j) + ".txt";
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
				cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

				// open position output file
				cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
				phiDisk = 0.5;
				// Initialze the system as disks
				cout << "	** Initializing at phiDisk = " << phiDisk << endl;
				cell_group.initializeGel(NV, phiDisk, sizedev, del);

				// change to DPM
				cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

				// set time scale
				cell_group.vertexDPMTimeScale(timeStepMag);



				// Compress then relax by FIRE
				cout << " Compress then relax by FIRE " << endl;

				//cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);
				double phiTargetTmp = 0.7;
				double deltaPhiTmp = 0.001;

				v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << NCELLS << "," << phiTargetTmp << endl;

				cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance);

				cellPacking2D jammed_state;
				cell_group.saveState(jammed_state);

				// output files
				extend = "_" + to_string(i) + to_string(j) + ".txt";
				//string positionF = "position" + extend;
				energyF = "energy" + extend;
				jammingF = "jam" + extend;
				lengthscaleF = "length" + extend;
				phiF = "phi" + extend;
				calAF = "calA" + extend;
				contactF = "contact" + extend;
				vF = "v" + extend;

				cell_group.loadState(jammed_state);
				cell_group.closeF();
				cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

				cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
				cell_group.activityCOM_brownian(T, v0, Dr, vtau, t_scale, frames);
				//cell_group.relaxF(Ktolerance, Ftolerance, Ptolerance);


			}
		}

		cout << "	** FINISHED **   " << endl;
	};


	void active_brownian()
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
		int NCELLS = 256;
		int NV = 16;
		int seed = 5;
		double Lini = 1.0;

		// activity
		double T = 1000.0;
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
		phiDisk = 0.5;
		cout << "	** Initializing at phiDisk = " << phiDisk << endl;
		cell_group.initializeGel(NV, phiDisk, sizedev, del);

		// change to DPM
		cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

		// set time scale
		cell_group.vertexDPMTimeScale(timeStepMag);



		// Compress then relax by FIRE
		cout << " Compress then relax by FIRE " << endl;

		double phiTargetTmp = 0.7;
		double deltaPhiTmp = 0.001;

		cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance);
		//cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);

		cellPacking2D jammed_state;
		cell_group.saveState(jammed_state);

		for (int i = 0; i < 1; i++) {
			for (int j = 0; j < 2; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;
				//v0 = 0.1;
				v0 = 1 + double(i) * 0.1;
				//Dr = 1.0 + double(j) * 1.0;
				Dr = 1e-3;
				kb = 0.0 + double(j) * 0.03;
				kl = 0.1;
				v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << NCELLS << endl;

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
				//cell_group.relaxF(Ktolerance, Ftolerance, Ptolerance);
				cell_group.activityCOM_brownian(T, v0, Dr, vtau, t_scale, frames);
				//cell_group.relaxF(Ktolerance, Ftolerance, Ptolerance);
			}
		}

		cout << "	** FINISHED **   " << endl;
	};

	void soft_particle_limit()
	{
		Ftolerance = 1e-7;
#pragma omp parallel
		{
			int ID = omp_get_thread_num();

			cout << "hello " << ID << endl;

		}

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");

#pragma omp parallel for
		for (int i = 0; i < 10; i++) {


			// system size
			int NCELLS = 16;
			int NV = 16;
			int seed = 5;
			double Lini = 1.0;

			// activity
			double T = 500.0;
			double v0;
			double Dr;
			double vtau = 1e-2;
			double t_scale = 1.00;
			double timeStepMag = 0.001;

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
			cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

			// open position output file
			cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
			double phiDisk = 0.45 + double(i) * 0.02;
			// Initialze the system as disks
			cout << "	** Initializing at phiDisk = " << phiDisk << endl;
			cell_group.initializeGel(NV, phiDisk, sizedev, del);

			// change to DPM
			cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

			// set time scale
			cell_group.vertexDPMTimeScale(timeStepMag);



			// Compress then relax by FIRE
			cout << " Compress then relax by FIRE " << endl;

			//cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);
			double phiTargetTmp = 0.7 + double(i) * 0.02;
			double deltaPhiTmp = 0.001;
			cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance);

			cellPacking2D jammed_state;
			cell_group.saveState(jammed_state);

			for (int j = 0; j < 10; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;

				//v0 = 0.1;
				v0 = 0.01 + double(j) * 0.01;
				//Dr = 1.0 + double(j) * 1.0;
				Dr = 1e-2;
				kl = 0.1;
				//kb = 0.0 + double(j) * 0.03;
				kb = 0.1;

				v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << NCELLS << "," << phiTargetTmp << endl;

				// output files
				extend = "_" + to_string(i) + to_string(j) + ".txt";
				//string positionF = "position" + extend;
				energyF = "energy" + extend;
				jammingF = "jam" + extend;
				lengthscaleF = "length" + extend;
				phiF = "phi" + extend;
				calAF = "calA" + extend;
				contactF = "contact" + extend;
				vF = "v" + extend;

				cell_group.loadState(jammed_state);
				cell_group.closeF();
				cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

				cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
				cell_group.activityCOM_brownian(T, v0, Dr, vtau, t_scale, frames);

			}
		}

		cout << "	** FINISHED **   " << endl;
	};

	void confluency()
	{	
#pragma omp parallel
		{
			int ID = omp_get_thread_num();

			cout << "hello " <<  ID << endl;

		}

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");

		// system size
		int NCELLS = 16;
		int NV = 16;
		int seed = 5;
		double Lini = 1.0;

		// activity
		double T = 10000.0;
		double Dr;
		double vtau = 1e-2;
		double t_scale = 1.00;

		//Dr = 1.0 + double(j) * 1.0;
		Dr = 1e-2;
		//kb = 0.0 + double(j) * 0.005;
		kb = 0.0;
		kl = 0.1;
		//kl = 0.05 + double(j) * 0.05;

		double phi_max = cal_phi_max(NCELLS, NV, seed, Lini, kl, kb);
		//double phi_max = 0.94;

#pragma omp parallel for
		for (int i = 0; i < 10; i++) {


			double calA0 = 1.12 + double(i) * 0.01;


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
			cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

			// open position output file
			cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
			phiDisk = 0.65;
			// Initialze the system as disks
			cout << "	** Initializing at phiDisk = " << phiDisk << endl;
			cell_group.initializeGel(NV, phiDisk, sizedev, del);

			// change to DPM
			cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

			// set time scale
			cell_group.vertexDPMTimeScale(timeStepMag);



			// Compress then relax by FIRE
			cout << " Compress then relax by FIRE " << endl;

			double phiTargetTmp = phi_max;
			double deltaPhiTmp = 0.001;

			cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance);
			//cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);

			cellPacking2D jammed_state;
			cell_group.saveState(jammed_state);

			for (int j = 0; j < 10; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;
	
				//v0 = 0.04;
				double v0 = 0.006 + double(j) * 0.003;

				v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << calA0 << "," << NCELLS << endl;

				// output files
				extend = "_" + to_string(i) + to_string(j) + ".txt";
				//string positionF = "position" + extend;
				energyF = "energy" + extend;
				jammingF = "jam" + extend;
				lengthscaleF = "length" + extend;
				phiF = "phi" + extend;
				calAF = "calA" + extend;
				contactF = "contact" + extend;
				vF = "v" + extend;

				cell_group.loadState(jammed_state);
				cell_group.closeF();
				cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

				cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
				cell_group.activityCOM_brownian(T, v0, Dr, vtau, t_scale, frames);

			}
		}

		cout << "	** FINISHED **   " << endl;
	};

	void confluency(char const *argv[])
	{	
#pragma omp parallel
		{
			int ID = omp_get_thread_num();

			cout << "hello " <<  ID << endl;

		}

		int index;
		string index_str = argv[1];
		stringstream indexss(index_str);
		indexss >> index;
		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");

		// system size
		int NCELLS = 16;
		int NV = 16;
		int seed = 5;
		double Lini = 1.0;

		// activity
		double T = 100000.0;
		int frames = 5000;
		double Dr;
		double vtau = 1e-2;
		double t_scale = 1.00;

		//Dr = 1.0 + double(j) * 1.0;
		Dr = 1e-2;
		//kb = 0.0 + double(j) * 0.005;
		kb = 0.0;
		kl = 0.1;
		//kl = 0.05 + double(j) * 0.05;

		//double phi_max = cal_phi_max(NCELLS, NV, seed, Lini, kl, kb);
		double phi_max = 0.95;

		for (int i = index; i < index + 1; i++) {


			double calA0 = 1.12 + double(i) * 0.01;


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
			cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

			// open position output file
			cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
			phiDisk = 0.65;
			// Initialze the system as disks
			cout << "	** Initializing at phiDisk = " << phiDisk << endl;
			cell_group.initializeGel(NV, phiDisk, sizedev, del);

			// change to DPM
			cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

			// set time scale
			cell_group.vertexDPMTimeScale(timeStepMag);



			// Compress then relax by FIRE
			cout << " Compress then relax by FIRE " << endl;

			double phiTargetTmp = phi_max;
			double deltaPhiTmp = 0.001;

			cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance);
			//cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);

			cellPacking2D jammed_state;
			cell_group.saveState(jammed_state);

#pragma omp parallel for
			for (int j = 0; j < 10; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;
	
				//v0 = 0.04;
				double v0 = 0.002 + double(j) * 0.002;

#pragma omp critical
				{
					v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << calA0 << "," << NCELLS << endl;
				}

				// output files
				string extend = "_" + to_string(i) + to_string(j) + ".txt";
				//string positionF = "position" + extend;
				string energyF = "energy" + extend;
				string jammingF = "jam" + extend;
				string lengthscaleF = "length" + extend;
				string phiF = "phi" + extend;
				string calAF = "calA" + extend;
				string contactF = "contact" + extend;
				string vF = "v" + extend;

				cellPacking2D local_cell_group;
				cell_group.saveState(local_cell_group);
				local_cell_group.closeF();
				local_cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

				local_cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
				local_cell_group.activityCOM_brownian(T, v0, Dr, vtau, t_scale, frames);

			}
		}

		cout << "	** FINISHED **   " << endl;
	};

	double cal_phi_max(int NCELLS, int NV, int seed, double Lini, double kl, double kb) {
		// output files
		string positionF = "position.txt";
		string energyF = "energy.txt";
		string jammingF = "jam.txt";
		string lengthscaleF = "length.txt";
		string phiF = "phi.txt";
		string calAF = "calA.txt";
		string contactF = "contact.txt";
		string vF = "v.txt";
		
		double calA0 = 1.18;

		// instantiate object
		cout << "	** Cell packing, NCELLS = " << NCELLS << endl;
		cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

		// open position output file
		cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

		// Initialze the system as disks
		phiDisk = 0.60;
		cout << "	** Initializing at phiDisk = " << phiDisk << endl;
		cell_group.initializeGel(NV, phiDisk, sizedev, del);

		// change to DPM
		cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

		// set time scale
		cell_group.vertexDPMTimeScale(timeStepMag);

		// Compress then relax by FIRE
		cout << " Compress then relax by FIRE " << endl;
		cell_group.findJamming(deltaPhi, Ftolerance, Ptolerance);

		return cell_group.packingFraction();

	};

	void const_ground_shape()
	{
#pragma omp parallel
		{
			int ID = omp_get_thread_num();

			cout << "hello " << ID << endl;

		}

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");

		// system size
		int NCELLS = 16;
		int NV = 16;
		int seed = 5;
		double Lini = 1.0;

		// activity
		double T = 10000.0;
		double Dr;
		double vtau = 1e-2;
		double t_scale = 1.00;

		//Dr = 1.0 + double(j) * 1.0;
		Dr = 1e-2;


		//double phi_max = cal_phi_max(NCELLS, NV, seed, Lini, kl, kb);
		double phi_max = 0.94;

		double ratio = 100.0;
		//ka = 10;

#pragma omp parallel for
		for (int i = 0; i < 10; i++) {


			double calA0 = 1.12;
			double kb = 0.0001 * (i + 1);
			double kl = ratio * kb;

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
			cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

			// open position output file
			cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
			phiDisk = 0.65;
			// Initialze the system as disks
			cout << "	** Initializing at phiDisk = " << phiDisk << endl;
			cell_group.initializeGel(NV, phiDisk, sizedev, del);

			// change to DPM
			cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

			// set time scale
			cell_group.vertexDPMTimeScale(timeStepMag);



			// Compress then relax by FIRE
			cout << " Compress then relax by FIRE " << endl;

			double phiTargetTmp = phi_max;
			double deltaPhiTmp = 0.001;

			cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance);
			//cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);

			cellPacking2D jammed_state;
			cell_group.saveState(jammed_state);

			for (int j = 0; j < 10; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;

				//v0 = 0.04;
				double v0 = 0.006 + double(j) * 0.003;

				v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << calA0 << "," << NCELLS << endl;

				// output files
				extend = "_" + to_string(i) + to_string(j) + ".txt";
				//string positionF = "position" + extend;
				energyF = "energy" + extend;
				jammingF = "jam" + extend;
				lengthscaleF = "length" + extend;
				phiF = "phi" + extend;
				calAF = "calA" + extend;
				contactF = "contact" + extend;
				vF = "v" + extend;

				cell_group.loadState(jammed_state);
				cell_group.closeF();
				cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

				cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
				cell_group.activityCOM_brownian(T, v0, Dr, vtau, t_scale, frames);

			}
		}

		cout << "	** FINISHED **   " << endl;
	};

	void const_ground_shape_arg(char const *argv[])
	{
#pragma omp parallel
		{
			int ID = omp_get_thread_num();

			cout << "hello " << ID << endl;

		}

		int index;
		string index_str = argv[1];
		stringstream indexss(index_str);
		indexss >> index;

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");

		// system size
		int NCELLS = 16;
		int NV = 16;
		int seed = 5;
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
			double calA0 = 1.18;
			//double kb = 0.0001 * (i + 1);
			//double kb = 0.00001 * pow(i + 1,2);
			double kb = 0.00005 * pow(i + 1,2);
			//double kl = ratio * kb;
			double kl = 0.1;

			double phi_max = 0.93;
			//double phi_max = 0.92;
			//double phi_max = cal_phi_max(NCELLS, NV, seed, Lini, kl, kb);

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
			cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

			// open position output file
			cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
			phiDisk = 0.65;
			//phiDisk = 0.60;
			// Initialze the system as disks
			cout << "	** Initializing at phiDisk = " << phiDisk << endl;
			cell_group.initializeGel(NV, phiDisk, sizedev, del);

			// change to DPM
			cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

			// set time scale
			cell_group.vertexDPMTimeScale(timeStepMag);



			// Compress then relax by FIRE
			cout << " Compress then relax by FIRE " << endl;

			double phiTargetTmp = phi_max;
			double deltaPhiTmp = 0.001;

			cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance);
			//cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);

			cellPacking2D jammed_state;
			cell_group.saveState(jammed_state);

#pragma omp parallel for 
			for (int j = 0; j < 10; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;

				//v0 = 0.04;
				//double v0 = 0.002 + double(j) * 0.002;
				double v0 = 0.0004 * double(i) + double(j+1) * 0.0015;
#pragma omp critical
				{
					v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << calA0 << "," << NCELLS << endl;
				}
				// output files
				string extend = "_" + to_string(i) + to_string(j) + ".txt";
				//string positionF = "position" + extend;
				string energyF = "energy" + extend;
				string jammingF = "jam" + extend;
				string lengthscaleF = "length" + extend;
				string phiF = "phi" + extend;
				string calAF = "calA" + extend;
				string contactF = "contact" + extend;
				string vF = "v" + extend;

				cellPacking2D local_cell_group;
				cell_group.saveState(local_cell_group);
				local_cell_group.closeF();
				local_cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

				local_cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
				local_cell_group.activityCOM_brownian(T, v0, Dr, vtau, t_scale, frames);

			}
		}

		cout << "	** FINISHED **   " << endl;
	};

	void cell_NVE_arg(char const *argv[])
	{
#pragma omp parallel
		{
			int ID = omp_get_thread_num();

			cout << "hello " << ID << endl;

		}

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
			//double kb = 0.0001 * (i + 1);
			//double kb = 0.00001 * pow(i + 1,2);
			double kb = 0.00001 * pow(i + 1,2);
			//double kl = ratio * kb;
			double kl = 0.1;

			double phi_max = 0.85;
			//double phi_max = 0.7;
			//double phi_max = cal_phi_max(NCELLS, NV, seed, Lini, kl, kb);

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
			cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

			// open position output file
			cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
			//phiDisk = 0.65;
			//phiDisk = 0.4;
			phiDisk = phi_max;
			// Initialze the system as disks
			cout << "	** Initializing at phiDisk = " << phiDisk << endl;
			cell_group.initializeGel(NV, phiDisk, sizedev, del);

			// change to DPM
			cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

			// set time scale
			cell_group.vertexDPMTimeScale(timeStepMag);



			// Compress then relax by FIRE
			cout << " Compress then relax by FIRE " << endl;

			double phiTargetTmp = phi_max;
			double deltaPhiTmp = 0.001;

			cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance);
			//cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);

			cellPacking2D jammed_state;
			cell_group.saveState(jammed_state);

#pragma omp parallel for 
			for (int j = 0; j < 10; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;

				//v0 = 0.04;
				//double v0 = 0.002 + double(j) * 0.002;
				//double v0 = 0.0004 * double(i) + double(j+1) * 0.0015;
				//double v0 = 0.0004 * double(i) + double(j+1) * 0.002;
				double v0 = double(j+1) * 0.0002;
#pragma omp critical
				{
					v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << calA0 << "," << NCELLS << endl;
				}
				// output files
				string extend = "_" + to_string(i) + to_string(j) + ".txt";
				//string positionF = "position" + extend;
				string energyF = "energy" + extend;
				string jammingF = "jam" + extend;
				string lengthscaleF = "length" + extend;
				string phiF = "phi" + extend;
				string calAF = "calA" + extend;
				string contactF = "contact" + extend;
				string vF = "v" + extend;

				cellPacking2D local_cell_group;
				cell_group.saveState(local_cell_group);
				local_cell_group.closeF();
				local_cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

				local_cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
				local_cell_group.cell_NVE(T, v0, Dr, vtau, t_scale, frames);

			}
		}

		cout << "	** FINISHED **   " << endl;
	};

	void random_NVE_arg(char const *argv[])
	{
#pragma omp parallel
		{
			int ID = omp_get_thread_num();

			cout << "hello " << ID << endl;

		}

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
			//double kb = 0.0001 * (i + 1);
			//double kb = 0.00001 * pow(i + 1,2);
			double kb = 0.00001 * pow(i + 1,2);
			//double kl = ratio * kb;
			double kl = 0.1;

			double phi_max = 0.9;
			//double phi_max = 0.4;
			//double phi_max = cal_phi_max(NCELLS, NV, seed, Lini, kl, kb);

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
			cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

			// open position output file
			cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
			//phiDisk = 0.65;
			//phiDisk = 0.4;
			phiDisk = phi_max;
			// Initialze the system as disks
			cout << "	** Initializing at phiDisk = " << phiDisk << endl;
			cell_group.initializeGel(NV, phiDisk, sizedev, del);

			// change to DPM
			cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

			// set time scale
			cell_group.vertexDPMTimeScale(timeStepMag);



			// Compress then relax by FIRE
			cout << " Compress then relax by FIRE " << endl;

			double phiTargetTmp = phi_max;
			double deltaPhiTmp = 0.001;

			cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance);
			//cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);

			cellPacking2D jammed_state;
			cell_group.saveState(jammed_state);

#pragma omp parallel for 
			for (int j = 0; j < 10; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;

				//v0 = 0.04;
				//double v0 = 0.002 + double(j) * 0.002;
				//double v0 = 0.0004 * double(i) + double(j+1) * 0.0015;
				//double v0 = 0.0004 * double(i) + double(j+1) * 0.002;
				double v0 = double(j+1) * 0.0002;
#pragma omp critical
				{
					v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << calA0 << "," << NCELLS << endl;
				}
				// output files
				string extend = "_" + to_string(i) + to_string(j) + ".txt";
				//string positionF = "position" + extend;
				string energyF = "energy" + extend;
				string jammingF = "jam" + extend;
				string lengthscaleF = "length" + extend;
				string phiF = "phi" + extend;
				string calAF = "calA" + extend;
				string contactF = "contact" + extend;
				string vF = "v" + extend;

				cellPacking2D local_cell_group;
				cell_group.saveState(local_cell_group);
				local_cell_group.closeF();
				local_cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

				local_cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
				local_cell_group.random_shape_NVE(T, v0, Dr, vtau, t_scale, frames);

			}
		}

		cout << "	** FINISHED **   " << endl;
	};

	void cell_NVE_probe_arg(char const *argv[])
	{
#pragma omp parallel
		{
			int ID = omp_get_thread_num();

			cout << "hello " << ID << endl;

		}

		int index;
		string index_str = argv[1];
		stringstream indexss(index_str);
		indexss >> index;

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");

		// system size
		int NCELLS = 16;
		int NV = 16;
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
			//double kb = 0.0001 * (i + 1);
			//double kb = 0.00001 * pow(i + 1,2);
			double kb = 0.00001 * pow(i + 1,2);
			//double kl = ratio * kb;
			double kl = 0.1;

			//double phi_max = 0.93;
			//double phi_max = 0.88;
			double phi_max = 0.7;
			//double phi_max = cal_phi_max(NCELLS, NV, seed, Lini, kl, kb);

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
			cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

			// open position output file
			cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
			//phiDisk = 0.65;
			//phiDisk = 0.63;
			phiDisk = 0.50;
			// Initialze the system as disks
			cout << "	** Initializing at phiDisk = " << phiDisk << endl;
			cell_group.initializeGel(NV, phiDisk, sizedev, del);

			// change to DPM
			cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

			// set time scale
			cell_group.vertexDPMTimeScale(timeStepMag);



			// Compress then relax by FIRE
			cout << " Compress then relax by FIRE " << endl;

			double phiTargetTmp = phi_max;
			double deltaPhiTmp = 0.001;

			cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance);
			//cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);

			cellPacking2D jammed_state;
			cell_group.saveState(jammed_state);

#pragma omp parallel for 
			for (int j = 0; j < 10; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;

				//v0 = 0.04;
				//double v0 = 0.002 + double(j) * 0.002;
				//double v0 = 0.0004 * double(i) + double(j+1) * 0.0015;
				//double v0 = 0.0004 * double(i) + double(j+1) * 0.002;
				double v0 = double(j+1) * 0.0002;
#pragma omp critical
				{
					v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << calA0 << "," << NCELLS << endl;
				}
				// output files
				string extend = "_" + to_string(i) + to_string(j) + ".txt";
				//string positionF = "position" + extend;
				string energyF = "energy" + extend;
				string jammingF = "jam" + extend;
				string lengthscaleF = "length" + extend;
				string phiF = "phi" + extend;
				string calAF = "calA" + extend;
				string contactF = "contact" + extend;
				string vF = "v" + extend;

				cellPacking2D local_cell_group;
				cell_group.saveState(local_cell_group);
				local_cell_group.closeF();
				local_cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

				local_cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
				local_cell_group.cell_NVE_probe(T, v0, Dr, vtau, t_scale, frames);

			}
		}

		cout << "	** FINISHED **   " << endl;
	};


	void sp_NVE_arg(char const* argv[])
	{
#pragma omp parallel
		{
			int ID = omp_get_thread_num();

			cout << "hello " << ID << endl;

		}

		int index;
		string index_str = argv[1];
		stringstream indexss(index_str);
		indexss >> index;

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");

		// system size
		int NCELLS = 16;
		int NV = 16;
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

		for (int i = index; i < index + 1; i++) {

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
			cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

			// open position output file
			cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
			//phiDisk = 0.7;
			phiDisk = 0.9;

			// Initialze the system as disks
			cout << "	** Initializing at phiDisk = " << phiDisk << endl;
			cell_group.initializeGel(NV, phiDisk, sizedev, del);

			// set time scale
			cell_group.vertexDPMTimeScale(timeStepMag);

			// Compress then relax by FIRE
			cout << " Compress then relax by FIRE " << endl;

			//double phiTargetTmp = phi_max;
			//double deltaPhiTmp = 0.001;

			//cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance);
			//cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);

			cellPacking2D jammed_state;
			cell_group.saveState(jammed_state);

#pragma omp parallel for 
			for (int j = 0; j < 10; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;

				//v0 = 0.04;
				//double v0 = 0.002 + double(j) * 0.002;
				//double v0 = 0.0004 * double(i) + double(j+1) * 0.0015;
				//double v0 = 0.0004 * double(i) + double(j + 1) * 0.002;
				double v0 = double(j+1) * 0.0002;
#pragma omp critical
				{
					v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << calA0 << "," << NCELLS << endl;
				}
				// output files
				string extend = "_" + to_string(i) + to_string(j) + ".txt";
				//string positionF = "position" + extend;
				string energyF = "energy" + extend;
				string jammingF = "jam" + extend;
				string lengthscaleF = "length" + extend;
				string phiF = "phi" + extend;
				string calAF = "calA" + extend;
				string contactF = "contact" + extend;
				string vF = "v" + extend;

				cellPacking2D local_cell_group;
				cell_group.saveState(local_cell_group);
				local_cell_group.closeF();
				local_cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
				local_cell_group.sp_NVE(T, v0, Dr, vtau, t_scale, frames);

			}
		}

		cout << "	** FINISHED **   " << endl;
	};

	void sp_NVE_probe_arg(char const* argv[])
	{
#pragma omp parallel
		{
			int ID = omp_get_thread_num();

			cout << "hello " << ID << endl;

		}

		int index;
		string index_str = argv[1];
		stringstream indexss(index_str);
		indexss >> index;

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");

		// system size
		int NCELLS = 64;
		int NV = 16;
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

		for (int i = index; i < index + 1; i++) {

			//double phi_max = 0.93;
			//double phi_max = 0.88;
			//double phi_max = 0.8;

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
			cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

			// open position output file
			cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
			phiDisk = 0.6;

			// Initialze the system as disks
			cout << "	** Initializing at phiDisk = " << phiDisk << endl;
			cell_group.initializeGel(NV, phiDisk, sizedev, del);

			// set time scale
			cell_group.vertexDPMTimeScale(timeStepMag);

			// Compress then relax by FIRE
			cout << " Compress then relax by FIRE " << endl;

			//double phiTargetTmp = phi_max;
			//double deltaPhiTmp = 0.001;

			//cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance);
			//cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);

			cellPacking2D jammed_state;
			cell_group.saveState(jammed_state);

#pragma omp parallel for 
			for (int j = 0; j < 10; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;

				//v0 = 0.04;
				//double v0 = 0.002 + double(j) * 0.002;
				//double v0 = 0.0004 * double(i) + double(j+1) * 0.0015;
				//double v0 = 0.0004 * double(i) + double(j + 1) * 0.002;
				double v0 = double(j+1) * 0.0002;
#pragma omp critical
				{
					v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << calA0 << "," << NCELLS << endl;
				}
				// output files
				string extend = "_" + to_string(i) + to_string(j) + ".txt";
				//string positionF = "position" + extend;
				string energyF = "energy" + extend;
				string jammingF = "jam" + extend;
				string lengthscaleF = "length" + extend;
				string phiF = "phi" + extend;
				string calAF = "calA" + extend;
				string contactF = "contact" + extend;
				string vF = "v" + extend;

				cellPacking2D local_cell_group;
				cell_group.saveState(local_cell_group);
				local_cell_group.closeF();
				local_cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
				local_cell_group.sp_NVE_probe(T, v0, Dr, vtau, t_scale, frames);

			}
		}

		cout << "	** FINISHED **   " << endl;
	};

	void bumpy_NVE_arg(char const* argv[])
	{
#pragma omp parallel
		{
			int ID = omp_get_thread_num();

			cout << "hello " << ID << endl;

		}

		int index;
		string index_str = argv[1];
		stringstream indexss(index_str);
		indexss >> index;

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");

		// system size
		int NCELLS = 16;
		int NV = 16;
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

		for (int i = index; i < index + 1; i++) {

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
			cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

			// open position output file
			cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
			//phiDisk = 0.7;
			phiDisk = 0.4;

			// Initialze the system as disks
			cout << "	** Initializing at phiDisk = " << phiDisk << endl;
			cell_group.initializeGel(NV, phiDisk, sizedev, del);

			// set time scale
			cell_group.vertexDPMTimeScale(timeStepMag);

			// Compress then relax by FIRE
			cout << " Compress then relax by FIRE " << endl;

			//double phiTargetTmp = phi_max;
			//double deltaPhiTmp = 0.001;

			//cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance);
			//cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);

#pragma omp parallel for 
			for (int j = 0; j < 10; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;

				//v0 = 0.04;
				//double v0 = 0.002 + double(j) * 0.002;
				//double v0 = 0.0004 * double(i) + double(j+1) * 0.0015;
				//double v0 = 0.0004 * double(i) + double(j + 1) * 0.002;
				double v0 = double(j + 1) * 0.0002;
#pragma omp critical
				{
					v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << calA0 << "," << NCELLS << endl;
				}
				// output files
				string extend = "_" + to_string(i) + to_string(j) + ".txt";
				//string positionF = "position" + extend;
				string energyF = "energy" + extend;
				string jammingF = "jam" + extend;
				string lengthscaleF = "length" + extend;
				string phiF = "phi" + extend;
				string calAF = "calA" + extend;
				string contactF = "contact" + extend;
				string vF = "v" + extend;

				cellPacking2D local_cell_group;
				cell_group.saveState(local_cell_group);
				local_cell_group.closeF();
				local_cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
				local_cell_group.bumpy_NVE(phiDisk, T, v0, Dr, vtau, t_scale, frames);

			}
		}

		cout << "	** FINISHED **   " << endl;
	};

	// main function
	void test()
	{
		string extend = "_jammed_" + to_string(0) + ".txt";
		//string positionF = "position" + extend;
		string energyF = "energy" + extend;
		string jammingF = "jam" + extend;
		string lengthscaleF = "length" + extend;
		string phiF = "phi" + extend;
		string calAF = "calA" + extend;
		string contactF = "contact" + extend;
		string vF = "v" + extend;

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");

		// system size
		int NCELLS = 2;
		int NV = 16;
		int seed = 5;
		double Lini = 1.0;

		// activity
		double T = 1000000.0;
		int frames = 50000;
		double v0 = 1.0;
		double Dr = 10.0;
		double vtau = 1e-2;
		double t_scale = 1.00;

		timeStepMag = 0.001;

		// instantiate object
		cout << "	** Cell packing, NCELLS = " << NCELLS << endl;
		cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

		// open position output file
		cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

		phiDisk = 0.2;
		// Initialze the system as disks
		cout << "	** Initializing at phiDisk = " << phiDisk << endl;
		cell_group.initializeGel(NV, phiDisk, sizedev, del);

		// change to DPM
		cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

		// set time scale
		cell_group.vertexDPMTimeScale(timeStepMag);



		// Compress then relax by FIRE
		cout << " Compress then relax by FIRE " << endl;

		double phiTargetTmp = 0.3;
		double deltaPhiTmp = 0.001;

		cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance);
		//cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);

		cellPacking2D jammed_state;
		cell_group.saveState(jammed_state);

		for (int i = 0; i < 1; i++) {
			for (int j = 0; j < 1; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;
				//v0 = 0.1;
				v0 = 0.0002+ double(i) * 0.0002;
				//Dr = 1.0 + double(j) * 1.0;
				Dr = 1e-2;
				//kl = 0.5;
				//ka = 10;
				//kb = 0.1 + double(j) * 0.01;
				//kb should be 0 ~ 0.03
				v0PrintObject << v0 << "," << Dr << "," << kb << endl;

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
				//cell_group.relaxF(Ktolerance, Ftolerance, Ptolerance);
				//cell_group.activityCOM_brownian(T, v0, Dr, vtau, t_scale, frames);
				//cell_group.relaxF(Ktolerance, Ftolerance, Ptolerance);
				cell_group.cell_NVE(T, v0, Dr, vtau, t_scale, frames);
			}
		}

		cout << "	** FINISHED **   " << endl;
	}


	void gravity()
	{
		int i = 0;
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

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");

		// system size
		int NCELLS = 2;
		int NV = 16;
		int seed = 5;
		double Lini = 1.0;

		gam = 0.1;
		kl = 0.1;
		const int NT = 2e5;			// number of time steps for flow simulation
		const int NPRINT = 1e3;			// number of steps between printing
		const double smallRadius 		= 0.5;			// radius fo smaller particles (diameter is length unit)
		const double sizeRatio 			= 1.4;			// ratio of small diameter to large diameter
		const double w0 				= 10.0;			// width of hopper reservoir (in units of small diameter)
		const double w 					= 1.5;			// orifice width (in units of small diameter)
		const double th 				= PI/4.0;		// hopper angle (pi - th = deflection angle from horizontal)
		const double phi0 				= 0.4;			// initial packing fraction
		double b = 0.1;

		Lini = 0.1*w0;

		timeStepMag = 0.01;

		vector<double> radii(NCELLS,0.0);
		for (int ci=0; ci<NCELLS; ci++){
			if (ci % 2 == 0)
				radii.at(ci) = smallRadius;
			else
				radii.at(ci) = smallRadius*sizeRatio;
		}

		// instantiate object
		cout << "	** Cell packing, NCELLS = " << NCELLS << endl;
		cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

		// open position output file
		cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

		phiDisk = 0.3;
		// Initialze the system as disks
		cout << "	** Initializing at phiDisk = " << phiDisk << endl;
		cell_group.initializeHopperDP(radii,w0,w,th,Lini,NV);

		// change to DPM
		cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

		// set time scale
		cell_group.vertexDPMTimeScale(timeStepMag);

		cellPacking2D jammed_state;
		cell_group.saveState(jammed_state);

		for (int i = 0; i < 1; i++) {
			for (int j = 0; j < 1; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;
				v0PrintObject << kl << "," << gam  << endl;
				
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
				//cell_group.relaxF(Ktolerance, Ftolerance, Ptolerance);
				double g = 0.1;
				cout << "	** Running hopper NVE with g = " << g << endl;
				// packingObject.hopperDPNVE(w0,w,th,g,T);
				cell_group.flowHopperDP(w0,w,th,g,b);
			}
		}

		cout << "	** FINISHED **   " << endl;
	};

	void Hopper(char const* argv[])
	{
#pragma omp parallel
		{
			int ID = omp_get_thread_num();

			cout << "hello " << ID << endl;

		}

		int index;
		string index_str = argv[1];
		stringstream indexss(index_str);
		indexss >> index;

		int i = index;
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

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");
		//std::ofstream clogPrintObject;
		//clogPrintObject.open("clog.txt");

		// system size
		int NCELLS = 64;
		int NV = 16;
		int seed = 5;
		double Lini = 1.0;

		//gam = 0.05;
		//kl = 0.1;
		//gam = 0.1;
		//kl = 0.5;
		kl = 5;
		ka = 10.0;
		const int NT = 1e6;			// number of time steps for flow simulation
		const int NPRINT = 1e3;			// number of steps between printing
		const double smallRadius = 0.5;			// radius fo smaller particles (diameter is length unit)
		//const double sizeRatio = 1.4;			// ratio of small diameter to large diameter
		const double sizeRatio = 1;
		const double w0 = 20.0;			// width of hopper reservoir (in units of small diameter)
		const double w = 1.5;			// orifice width (in units of small diameter)
		const double th = PI / 4.0;		// hopper angle (pi - th = deflection angle from horizontal)
		const double phi0 = 0.4;			// initial packing fraction
		double b = 0.1;

		Lini = 0.1 * w0;

		//timeStepMag = 0.01;

		vector<double> radii(NCELLS, 0.0);
		for (int ci = 0; ci < NCELLS; ci++) {
			if (ci % 2 == 0)
				radii.at(ci) = smallRadius;
			else
				radii.at(ci) = smallRadius * sizeRatio;
		}

		// instantiate object
		cout << "	** Cell packing, NCELLS = " << NCELLS << endl;
		seed = i;
		cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

		// open position output file
		cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

		cell_group.initializeHopperDP(radii, w0, w, th, Lini, NV);

		// change to DPM
		cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

		// set time scale
		cell_group.vertexDPMTimeScale(timeStepMag);

		#pragma omp parallel for 
			for (int j = 0; j < 10; j++) {
				
				cout << "Loop i, j = " << i << "," << j << endl;
				//double gam = 0.05 + 0.02 * (j);
				double gam = 0.5 + 0.2 * (j);
				double g = 0.05;
				
				cout << "	** Running hopper NVE with g = " << g << endl;
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

				cellPacking2D local_cell_group;
				cell_group.saveState(local_cell_group);
				//local_cell_group.closeF();
				local_cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

				local_cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
				//cell_group.relaxF(Ktolerance, Ftolerance, Ptolerance);
				
				
				// packingObject.hopperDPNVE(w0,w,th,g,T);
				int result = local_cell_group.flowHopperDP(w0, w, th, g, b);
		#pragma omp critical
			{	
				v0PrintObject << kl << "," << gam << "," << g << "," << result << endl;
			}			
			}
		cout << "	** FINISHED **   " << endl;
		//clogPrintObject.close();
		v0PrintObject.close();
	};

	void  Hopper_width(char const* argv[])
	{
		int index_i, index_j;
		string index_i_str = argv[1];
		string index_j_str = argv[2];
		stringstream index_i_ss(index_i_str);
		stringstream index_j_ss(index_j_str);
		index_i_ss >> index_i;
		index_j_ss >> index_j;

		cout << "Loop i, j = " << index_i << "," << index_j << endl;

		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");
		//std::ofstream clogPrintObject;
		//clogPrintObject.open("clog.txt");

		// system size
		int NCELLS = 64;
		int NV = 16;
	
		double Lini = 1.0;

		//gam = 0.05;
		//kl = 0.1;
		kl = 1.0;
		ka = 1.0;
		const int NT = 5e6;			// number of time steps for flow simulation
		const int NPRINT = 1e3;			// number of steps between printing
		const double smallRadius = 0.5;			// radius fo smaller particles (diameter is length unit)
		const double sizeRatio = 1.4;			// ratio of small diameter to large diameter
		//const double sizeRatio = 1;
		const double w0 = 20.0;			// width of hopper reservoir (in units of small diameter)
		
		const double th = PI / 4.0;		// hopper angle (pi - th = deflection angle from horizontal)
		const double phi0 = 0.4;			// initial packing fraction
		double b = 0.1;

		Lini = 0.1 * w0;

		//timeStepMag = 0.01;

		vector<double> radii(NCELLS, 0.0);
		for (int ci = 0; ci < NCELLS; ci++) {
			if (ci % 2 == 0)
				radii.at(ci) = smallRadius;
			else
				radii.at(ci) = smallRadius * sizeRatio;
		}

		double w_scale = 0.5 + 0.1 * index_j;			// orifice width (in units of mean diameter)
		//double w_scale = 0.5 + 0.05 * index_j;
		//double gam = 0.5 + 0.2 * (j);
		double gam = 0;
		//double g = 0.05;
		double kl = 5, g = 0.02;

		double w = w_scale * (1 + sizeRatio) / 2;
		// output files
		string extend = "_jammed_" + to_string(index_i) + ".txt";
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
		int seed = index_i;
		cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

		// open position output file
		cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

		cell_group.initializeHopperDP(radii, w0, w, th, Lini, NV);

		// change to DPM
		cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

		// set time scale
		cell_group.vertexDPMTimeScale(timeStepMag);
		
		cout << "	** Running hopper NVE with g = " << g << endl;
		extend = "_" + to_string(index_i) + to_string(index_j) + ".txt";
		
		// output files
		//string positionF = "position" + extend;
		energyF = "energy" + extend;
		jammingF = "jam" + extend;
		lengthscaleF = "length" + extend;
		phiF = "phi" + extend;
		calAF = "calA" + extend;
		contactF = "contact" + extend;
		vF = "v" + extend;

		
		cell_group.closeF();
		cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

		//cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
		//cell_group.relaxF(Ktolerance, Ftolerance, Ptolerance);
		
		
		// packingObject.hopperDPNVE(w0,w,th,g,T);
		int result = cell_group.flowHopperDP(w0, w, th, g, b);

		v0PrintObject << kl << "," << gam << "," << g << "," << w_scale << "," << result << endl;
		
		cout << "	** FINISHED **   " << endl;
		//clogPrintObject.close();
		v0PrintObject.close();
	};


	void test_ground_shape()
	{
		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");
		int frames = 10;

		// system size
		int NCELLS = 2;
		int NV = 16;
		int seed = 5;
		double Lini = 1.0;

		// activity
		double T = 1.0;
		double Dr;
		double vtau = 1e-2;
		double t_scale = 1.00;

		//Dr = 1.0 + double(j) * 1.0;
		Dr = 1e-2;
		//kb = 0.0 + double(j) * 0.005;

		//kl = 0.05 + double(j) * 0.05;

		double phi_max = 0.3;


		for (int i = 0; i < 10; i++) {
			ka = 10;
			double calA0 = 1.12 ;
			kl = 0.05 + double(i) * 0.05;
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
			cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

			// open position output file
			cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
			phiDisk = 0.3;
			// Initialze the system as disks
			cout << "	** Initializing at phiDisk = " << phiDisk << endl;
			cell_group.initializeGel(NV, phiDisk, sizedev, del);

			// change to DPM
			cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

			// set time scale
			cell_group.vertexDPMTimeScale(timeStepMag);



			// Compress then relax by FIRE
			cout << " Compress then relax by FIRE " << endl;

			//cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);
			double phiTargetTmp = phi_max;
			double deltaPhiTmp = 0.001;
			cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance);

			cellPacking2D jammed_state;
			cell_group.saveState(jammed_state);

			for (int j = 0; j < 10; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;

				//v0 = 0.04;
				double v0 = 0.001;
				kb = 0.0 + double(j) * 0.0002;
				v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << calA0 << "," << NCELLS << endl;

				// output files
				extend = "_" + to_string(i) + to_string(j) + ".txt";
				//string positionF = "position" + extend;
				energyF = "energy" + extend;
				jammingF = "jam" + extend;
				lengthscaleF = "length" + extend;
				phiF = "phi" + extend;
				calAF = "calA" + extend;
				contactF = "contact" + extend;
				vF = "v" + extend;

				cell_group.loadState(jammed_state);
				cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

				double Fcheck, Kcheck;
				cell_group.fireMinimizeF(Ftolerance, Fcheck, Kcheck);

				cell_group.closeF();
				cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

				
				cell_group.activityCOM_brownian(T, v0, Dr, vtau, t_scale, frames);

			}
		}

		cout << "	** FINISHED **   " << endl;
	};

	void mesaure_ground_shape()
	{
		std::ofstream v0PrintObject;
		v0PrintObject.open("v0.txt");
		int frames = 10;

		// system size
		int NCELLS = 2;
		int NV = 16;
		int seed = 5;
		double Lini = 1.0;

		// activity
		double T = 1.0;
		double Dr;
		double vtau = 1e-2;
		double t_scale = 1.00;

		//Dr = 1.0 + double(j) * 1.0;
		Dr = 1e-2;
		//kb = 0.0 + double(j) * 0.005;

		//kl = 0.05 + double(j) * 0.05;

		double phi_max = 0.3;


		for (int i = 0; i < 10; i++) {
			ka = 1;
			double calA0 = 1.12 ;
			double ratio = 100;
			double kb = 0.0001 * (i + 1);
			double kl = ratio * kb;
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
			cellPacking2D cell_group(NCELLS, NT, NPRINT, Lini, seed);

			// open position output file
			cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);
			phiDisk = 0.3;
			// Initialze the system as disks
			cout << "	** Initializing at phiDisk = " << phiDisk << endl;
			cell_group.initializeGel(NV, phiDisk, sizedev, del);

			// change to DPM
			cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

			// set time scale
			cell_group.vertexDPMTimeScale(timeStepMag);



			// Compress then relax by FIRE
			cout << " Compress then relax by FIRE " << endl;

			//cell_group.findJamming(deltaPhi, Ktolerance, Ftolerance, Ptolerance);
			double phiTargetTmp = phi_max;
			double deltaPhiTmp = 0.001;
			cell_group.qsIsoCompression(phiTargetTmp, deltaPhiTmp, Ftolerance);

			cellPacking2D jammed_state;
			cell_group.saveState(jammed_state);

			for (int j = 0; j < 1; j++) {

				cout << "Loop i, j = " << i << "," << j << endl;

				//v0 = 0.04;
				double v0 = 0.0001;
				v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << calA0 << "," << NCELLS << endl;

				// output files
				extend = "_" + to_string(i) + to_string(j) + ".txt";
				//string positionF = "position" + extend;
				energyF = "energy" + extend;
				jammingF = "jam" + extend;
				lengthscaleF = "length" + extend;
				phiF = "phi" + extend;
				calAF = "calA" + extend;
				contactF = "contact" + extend;
				vF = "v" + extend;

				cell_group.loadState(jammed_state);
				cell_group.forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);

				double Fcheck, Kcheck;
				cell_group.fireMinimizeF(Ftolerance, Fcheck, Kcheck);

				cell_group.closeF();
				cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF);

				
				cell_group.activityCOM_brownian(T, v0, Dr, vtau, t_scale, frames);

			}
		}

		cout << "	** FINISHED **   " << endl;
	};








































};


#endif
