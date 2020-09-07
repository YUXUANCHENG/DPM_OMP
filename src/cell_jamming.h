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
			double kl = ratio * kb;

			double phi_max = 0.94;
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

	void test()
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
		int NCELLS = 2;
		int NV = 16;
		int seed = 5;
		double Lini = 1.0;

		// activity
		double T = 10.0;
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
				v0 = 0.1 + double(i) * 0.1;
				//Dr = 1.0 + double(j) * 1.0;
				Dr = 1e-2;
				kl = 0.5;
				//ka = 10;
				kb = 0.1 + double(j) * 0.01;
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
				cell_group.activityCOM_brownian(T, v0, Dr, vtau, t_scale, frames);
				//cell_group.relaxF(Ktolerance, Ftolerance, Ptolerance);
			}
		}

		cout << "	** FINISHED **   " << endl;
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
