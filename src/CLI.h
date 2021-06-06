#ifndef CLI_H
#define CLI_H

// include files
#include "cellPacking2D.h"
#include "bumpy.h"
#include "deformableParticles2D.h"
#include "bumpyEllipse.h"
#include "bumpyDimer.h"
#include <sstream>
#include <cmath>
#include <omp.h>

template <class Ptype = cellPacking2D>
class DPM_CLI {
public:
	// define PI
	const double PI = 4.0 * atan(1);

	// length paramaters
	int NT = 1e7;
	int NPRINT = 10000;

	// simulation constants
	double sizedev = 0.1;			        // std dev of cell sizes
	double timeStepMag = 0.005;		// time step in MD units (zeta * lenscale / forcescale)

	// disk constants
	double phiDisk = 0.75;			// initial packing fraction of disks

	// compression constants
	// const double phiTarget = 1.03;			// cell packing fraction (regardless of final pressure)
	const double deltaPhi = 0.001;		// compression step size

	// force parameters
	double kl = 0.1;				// perimeter force constant
	double ka = 1.0;				// area force constant
	double kb = 0.0;				// bending energy constant
	double gam = 0.0;				// surface tension force constant

	double kint = 1.0;				// interaction energy constant
	const double del = 1.0;				// width of vertices in units of l0, vertex sep on regular polygon
	const double aInitial = 0.0;				// attraction parameter to start

	// ratio of preferred perimeter^2 to preferred area
	double calA0 = 1.2;

	// tolerances
	double Ftolerance = 1e-10;			// force tolerance (for FIRE min)
	double Ptolerance = 1e-8;
	double Ktolerance = 1e-16;
	string extend;

	// system size
	int NCELLS = 64;
	int NV = 16;
	int seed = 1;
	double Lini = 1.0;
	int NBx = 8;
	int NBy = NBx;

	double Phi_to_PhiJ = 0.03;

	// activity
	double T = 10000.0;
	int frames = 20000;
	double Dr;
	double vtau = 1e-2;
	double t_scale = 1.00;
	std::ofstream v0PrintObject;
	int index_i, index_j;

	cellPacking2D* particles;

	void _createParticles(char const* argv[])
	{
		
		setIndex(argv);
		//seed = index_i;

		v0PrintObject.open("v0.txt");
		setKB();
		setSeed();
		setPhiDisk();

		// output files
		string extend = "_jammed_" + to_string(index_i) + ".txt";
		string energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF;
		produceFileName(extend, energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF);
		// instantiate object
		cout << "	** Cell packing, NCELLS = " << NCELLS << endl;
		particles = new Ptype(NCELLS, NT, NPRINT, Lini, seed);
		std::cout << typeid(particles).name() << '\n';
		particles->openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF);

	}
	virtual void setPhiDisk(){
		//phiDisk = 0.75;
		phiDisk = 0.70 + index_i * 0.02;
	}
	virtual void setSeed() {
		//seed = index_i;
		seed = 0;
	}
	virtual void setKB() {
		//double ratio = 100.0;
		kb = 0.00001 * pow(index_i + 1, 2);
		//kb = 0.0;
		//phiDisk = 0.8 + 0.01 * index_i;
		//double kl = ratio * kb;
	}

	virtual void setIndex(char const* argv[]) {
		string index_str = argv[1];
		stringstream indexss(index_str);
		indexss >> index_i;
	}

	void produceFileName(string extend, string& energyF, string& jammingF, string& lengthscaleF,
		string& phiF, string& calAF, string& contactF, string& vF,  string& ISF) {
		//string positionF = "position" + extend;
		energyF = "energy" + extend;
		jammingF = "jam" + extend;
		lengthscaleF = "length" + extend;
		phiF = "phi" + extend;
		calAF = "calA" + extend;
		contactF = "contact" + extend;
		vF = "v" + extend;
		ISF = "isf" + extend;
	}

	virtual void prepareSystem() {
		// Initialze the system as disks
		cout << "	** Initializing at phiDisk = " << phiDisk << endl;
		particles->initializeGel(NV, phiDisk, sizedev, del);
		particles->initialize_subsystems(NBx, NBy);
		particles->forceVals(calA0, kl, ka, gam, kb, kint, del, aInitial);
		particles->vertexDPMTimeScale(timeStepMag);
		particles->compressToInitial(phiDisk, deltaPhi, Ftolerance);
	}

	void findJamming(char const* argv[]) {
		qscompress(argv);
		particles->findJamming(deltaPhi, Ftolerance, Ptolerance);
	}

	void qscompress(char const* argv[]) {
		createParticles(argv);
		prepareSystem();
	}

	void toDeltaPhi(double delta) {
		double phi = particles->getphi();
		particles->qsIsoCompression(phi + delta, deltaPhi, Ftolerance);
	}

	void _NVE() {
//#pragma omp parallel for 
		for (int j = 0; j < 10; j++) {

			cout << "Loop i, j = " << index_i << "," << j << endl;
			//double v0 = 0.0004 * double(i) + double(j+1) * 0.0015;
			//double v0 = 0.0004 * double(i) + double(j + 1) * 0.002;
			double v0 = double(j + 1) * 0.0002;
#pragma omp critical
			{
				v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << calA0 << "," << NCELLS << endl;
			}
			// output files
			string extend = "_" + to_string(index_i) + to_string(j) + ".txt";
			string energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF;
			produceFileName(extend, energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF);

			Ptype local_cell_group;
			particles->saveState(local_cell_group);
			local_cell_group.closeF();
			local_cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF);
			//local_cell_group.NVEsimulation(T, v0, t_scale, frames);
			local_cell_group.initialize_subsystems();
			local_cell_group.NVEsimulation(T, v0, t_scale, frames);
			//local_cell_group.LangevinSimulation(T, v0, t_scale, frames);
		}
	}

	virtual void createParticles(char const* argv[])
	{
		_createParticles(argv);
	}

	virtual void NVEvsDPhi(char const* argv[])
	{
		findJamming(argv);
		toDeltaPhi(Phi_to_PhiJ);
		_NVE();
	}

	virtual void NVE(char const* argv[])
	{
		qscompress(argv);
		_NVE();
	}

	virtual void calTao(char const* argv[])
	{
		qscompress(argv);

		double preset_time = T;

		for (int j = 0; j < 10; j++) {
			cout << "Loop i, j = " << index_i << "," << j << endl;

			// double start = -8.5 + 0.25 * index_i;
			// double interval = 0.25 - 0.02 * index_i;
			double start = -7.8 + 0.25 * index_i;
			double interval = 0.25 - 0.02 * index_i;
			double v0 = exp(start + interval * (9 - j));

			// output files
			string extend = "_" + to_string(index_i) + to_string(j) + ".txt";
			string energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF;
			produceFileName(extend, energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF);

			Ptype local_cell_group;
			particles->saveState(local_cell_group);
			local_cell_group.closeF();
			local_cell_group.openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF);
			//local_cell_group.sp_NVE(T, v0, Dr, vtau, t_scale, frames);
			local_cell_group.initialize_subsystems();
			if (j == 0) local_cell_group.LangevinSimulation(4000, v0, t_scale, 10);
			double* result = local_cell_group.NVE_tao(preset_time, v0, Dr, vtau, t_scale, frames);
			v0PrintObject << v0 << "," << Dr << "," << kb << "," << kl << "," << calA0 << "," << NCELLS << "," << result[0] << "," << result[1] << endl;
			preset_time = result[0] * 10;
			delete[] result;
			local_cell_group.saveState(*particles);
		}
	}
};

template <class Ptype = Bumpy>
class Bumpy_CLI : public DPM_CLI<Ptype> {
public:
	//typedef Bumpy particleType;

	virtual void prepareSystem() {
		// Initialze the system as disks
		cout << "	** Initializing at phiDisk = " << this->phiDisk << endl;
		this->particles->initializeGel(this->NV, this->phiDisk, this->sizedev, this->del);
		this->particles->initialize_subsystems(this->NBx, this->NBy);
		this->particles->vertexDPMTimeScale(this->timeStepMag);
		this->particles->compressToInitial(this->phiDisk, this->deltaPhi, this->Ftolerance);
	}

};

template <class Ptype = BumpyEllipse>
class BumpyEllipse_CLI : public Bumpy_CLI<Ptype> {
public:
	//typedef BumpyEllipse particleType;
	double ratio = 1.7;

	virtual void setPhiDisk(){
		this->phiDisk = 0.86;
	}

	virtual void prepareSystem() {
		// Initialze the system as disks
		cout << "	** Initializing at phiDisk = " << this->phiDisk << endl;
		this->particles->setRatio(this->ratio);
		this->particles->initializeGel(this->NV, this->phiDisk, this->sizedev, this->del);
		this->particles->vertexDPMTimeScale(this->timeStepMag);
		this->particles->compressToInitial(this->phiDisk, this->deltaPhi, this->Ftolerance);
	}

};

template <class Ptype = BumpyDimer>
class BumpyDimer_CLI : public BumpyEllipse_CLI<Ptype> {
public:
	//typedef BumpyDimer particleType;

	virtual void setPhiDisk(){
		//phiDisk = 0.86;
		this->seed = 1;
		this->phiDisk = 0.65 + 0.01 * this->index_i;
		this->ratio = 1.6;
		this->NV = 16;
	}

};

#endif