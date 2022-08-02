#ifndef CELLPACKING2D_H
#define CELLPACKING2D_H


/*

	-- PACKING OF CELLS CLASS -- 

	* Class for systems of deformable particles, 
		holds deformable particles in packing.

	* Mostly wrapper functions for objects in cells array

	* Information on temperature, boundary changes (extension, compression, shear),

	* Note: everything in d=2, not generalizable to any d

*/
#include "subspace.h"
#include "deformableParticles2D.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <random>
#include <map>


// constants
const double PI = 4 * atan(1);

//class subspace;
//class cvpair;

class cellPacking2D{
public:

	frictionlessSubspace* subsystem = nullptr;
	vector<double> diskRadii;
	// std::map<cvpair, std::vector<cvpair>> collisionMap;
	std::vector<int> N_systems;
	virtual void initialize_subsystems() { ; };
	virtual void initialize_subsystems(int N_x, int N_y){ ; };
	virtual VECTOR6 gradient(vector<VECTOR2> & v, vector<VECTOR2> & e, double eps, int mode){return VECTOR6();};
	

	// int scalars
	int NDIM;						// spatial dimension (will always be 2)
	int NCELLS;						// number of cells
	int NT;							// number of total time steps
	int NPRINT;						// number of time steps between print steps

	// contact numbers
	int Ncc;						// number of cell-cell contacts
	int Nvv;						// number of vertex-vertex contacts

	// double scalars
	double seed;					// initial seed
	double dt;						// time step size
	double dt0;						// initial time step size (for dynamic time step methods)
	double T;						// system temperature
	double phi;						// system packing/volume fraction
	double shearStrain;				// applied shear strain to compute shear modulus
	int gDire = 0;
	int gOn = 1;
	double cutoff = 0.2;
	double insideCutoffFactor = 1;

	// boundary lengths
	std::vector<double> L;
	std::vector<double> BoundaryCoor;
	double pistonX, pistonY;

	// virial stresses
	double sigmaXX, sigmaXY, sigmaYX, sigmaYY;
	int print_frequency = 1;

	
	// array of cells
	deformableParticles2D* cellArray;
	double coefu = 0;
	double coefv = 0;
	double spK = 10;
	

	// contact matrix
	int* contactMatrix;				// array of contacts between particles (= number of vv contacts between two cells)

	// ofstream objects
	std::ofstream packingPrintObject;
	std::ofstream energyPrintObject;
	std::ofstream statPrintObject;
	std::ofstream jamPrintObject;
	std::ofstream lengthscalePrintObject;
	std::ofstream phiPrintObject;
	std::ofstream calAPrintObject;
	std::ofstream contactPrintObject;
	std::ofstream vPrintObject;
	std::ofstream ISFPrintObject;


	// Constructors and Destructors
	void defaultvars();
	cellPacking2D();
	cellPacking2D(int ncells, int nt, int nprint, double l, double s);
	cellPacking2D(int ncells, int ntumor, int tumorNV, int adiposeNV, double tumorCalA, double adiposeCalA, int s);
	cellPacking2D(std::string& inputFile, double T0, double s);
	~cellPacking2D();

	// new interface

	virtual void printRoutine(int count, int print_frequency, double t, double init_E, double init_U);
	virtual void printSubRoutine(int count, int print_frequency);
	virtual void NVEsimulation(double T, double v0, double t_scale, int frames);
	virtual void LangevinSimulation(double T, double v0, double t_scale, int frames);
	virtual int hopperSimulation(double w0, double w, double th, double g, double b);
	void pistonForce(double w0, double w, double th, double pressure);
	int forceWithBead(int ci, int vi, std::vector<double> bead);
	virtual void resetV();
	virtual void setRatio(double ratio) { ; };
	void calInertia();
	double calOriginalHeight();
	double calHeight();
	void changeL0(double factor1, double factor2);
	double calContactAng();
	double calContactLength();
	double matchPreset(double presetCalA, double presetAngle, double threashold);
	//bool matchPreset(double presetCalA, double presetAngle, double threashold);



	// operators
	void operator=(cellPacking2D& onTheRight);	// assign one configuration to another object
	void saveState(cellPacking2D& saveObject);
	void loadState(cellPacking2D& saveObject);

	// initialize sizes
	void initializeBidisperse(int NV, double phi0, double sizeRatio, double sizeFraction, double delval);

	// initialize velocities
	void initializeVelocities(double tmp0);

	// file openers
	void openPackingObject(std::string& str){
		packingPrintObject.open(str.c_str());
		if (!packingPrintObject.is_open()) {
			std::cout << "	ERROR: packingPrintObject could not open " << str << "..." << std::endl;
			exit(1);
		}
	}

	void openEnergyObject(std::string& str){
		energyPrintObject.open(str.c_str());
		if (!energyPrintObject.is_open()) {
			std::cout << "	ERROR: energyPrintObject could not open " << str << "..." << std::endl;
			exit(1);
		}
	}

	void openStatObject(std::string& str){
		statPrintObject.open(str.c_str());
		if (!statPrintObject.is_open()) {
			std::cout << "	ERROR: statPrintObject could not open " << str << "..." << std::endl;
			exit(1);
		}
	}
	
	void openJamObject(std::string& str, std::string& str1, std::string& str2, std::string& str3, std::string& str4, std::string& str5, std::string& str6) {
		jamPrintObject.open(str.c_str());
		lengthscalePrintObject.open(str1.c_str());
		phiPrintObject.open(str2.c_str());
		calAPrintObject.open(str3.c_str());
		contactPrintObject.open(str4.c_str());
		vPrintObject.open(str5.c_str());
		ISFPrintObject.open(str6.c_str());
		if (!jamPrintObject.is_open()) {
			std::cout << "	ERROR: jamPrintObject could not open " << str << "..." << std::endl;
			exit(1);
		}
	}

	// read in initial configuration from file
	void readInFromFile(std::ifstream& inputFileObject);

	// getters
	int getNCELLS() { return NCELLS; };
	int getNT() { return NT; };
	int getNPRINT() { return NPRINT; };
	double getdt() { return dt; };
	double getT() { return T; };
	double getphi() { return phi; };
	double getShearStrain() { return shearStrain; };
	double getSigmaXX() { return sigmaXX; };
	double getSigmaXY() { return sigmaXY; };
	double getSigmaYX() { return sigmaYX; };
	double getSigmaYY() { return sigmaYY; };
	void setFriction(double coefu, double coefv){
		this->coefu = coefu;
		this->coefv = coefv;
	};
	void setSPk(double spK){
		this->spK = spK;
	};

	// box len
	double getL(int d) { return L.at(d); };

	// array of cell objects
	deformableParticles2D& cell(int ci);		// return cell object ci
	int nframes();								// number of frames in the simulation

	// contact information
	int cmindex(int ci, int cj);				// contact matrix index
	int contacts(int ci, int cj);				// contact matrix element
	int vvContacts();							// total number of vertex-vertex contacts
	int vvContacts(int ci);						// total number of vertex-vertex contacts on cell ci
	int ccContacts();							// total number of cell-cell contacts
	int ccContacts(int ci);						// total number of cell-cell contacts on cell ci
	void setContact(int ci, int cj, int val);	// set contact to a certain value
	void addContact(int ci, int cj);			// add a single contact
	void deleteContact(int ci, int cj);			// delete a single contact
	void resetContacts();						// reset all contacts to 0

	// system-wide calculations
	int totalNumberOfContacts();			// calculate total number of contacts
	double packingFraction();				// calculate system packing fraction
	double totalPotentialEnergy();
	double totalPotentialEnergy_probe();				// calculate total potential energy
	double interactionPotentialEnergy();	// calculate only interaction potential energy
	double totalKineticEnergy();
	double totalKineticEnergy_probe();			// calculate total kinetic energy
	double maxForceMagnitude();				// calculate magnitude of largest force
	double maxNetForceMagnitude();			// calculate magnitude of largest net force on any cell
	double forceRMS(); 						// RMS net force
	double meanAsphericity();				// calculate system-average asphericity

	// setters
	void setNCELLS(int nc) { NCELLS = nc; };
	void setNT(int nt) { NT = nt; };
	void setNPRINT(int nprint) { NPRINT = nprint; };
	void setT(double val) { T = val; };
	void setL(int d, double val) { L.at(d) = val; };
	void setdt(double val) { dt0 = val; dt = dt0; };
	void setShearStrain(double val) { shearStrain = val; };
	void updatePackingFraction() { phi = packingFraction(); };

	// set force values for all cells to be the same
	virtual void forceVals(double calA0, double kl, double ka, double gam, double kb, double kint, double del, double a);

	void vertexDPMTimeScale(double timeStepMag);
	void setPackingFraction(double val);
	void setAsphericity(double val);
	void setAsphericity(int ci, double val);
	virtual void scaleLengths(double val);
	void rescaleVelocities(double temperature);
	void rescaleVelocities_probe(double temperature);
	int removeRattlers(int krcrs);

	/**************************

		Forces and position 
			updates

	***************************/

	virtual void calculateForces();
	void gelationForces();


	/**************************

		FIRE energy minimzation

	***************************/

	// FIRE 2.0 relaxation functions
	void fireMinimizeP(double Ptol, double Ktol);
	virtual void fireMinimizeF(double Ftol, double& Ftest, double& Ktest);


	/**************************

		Simulation Functions:

	***************************/

	// NVE test function
	void cellNVE();
	void cellOverDamped();

	// Find jammed states
	void findJamming(double dphi0, double Ftol, double Ptol);
	void enthalpyMin(double dphi0, double Ftol, double Ptol);

	// compress isotropically to a target packing fraction
	virtual void qsIsoCompression(double phiTarget, double deltaPhi, double Ftol);
	virtual void compressToInitial(double phiTarget, double deltaPhi, double Ftol);

	// compute the instantaneous shear modulus
	double shearModulus();

	// compute the vibrational density of states
	void vdos();

	// Gelation functions
	void twoParticleContact(int NV);
	virtual void initializeGel(int NV, double phiDisk, double sizeDispersion, double delval);
	void qsIsoGelRatchet(double phiGel, double deltaPhi, double plThresh, double dl0, double calA0max, double timeStepMag);
	void ratchetPerimeter(double plThresh, double dl0, double calA0max);
	void cellRK4();







	// Sticky SP particle functions
	void initializeStickySP(std::vector<double>& radii, double phiDisk, double sizeDispersion);
	void stickySPTriangularLattice(std::vector<double>& radii, double phiDisk);
	void stickySPGelationQS(std::vector<double>& radii, double phiGel, double dphiGel, double attractiveParam);
	void stickySPGelationRate(std::vector<double>& radii, double phiGel, double gelRate, double attractiveParam, double timeStepMag);
	void printPositionsStickySP(std::vector<double>& radii);
	void printEnergyStickySP();


	// Repulsive SP particle functions
	void fireMinimizeSP(std::vector<double>& lenscales);
	void fireMinimizeSP(std::vector<double>& radii, double attractiveParam);
	void fireMinimize_disk(vector<double>& lenscales);
	void fireMinimize_bummpy();
	void spForces(std::vector<double>& lenscales);
	void spAttractiveForces(std::vector<double>& radii, double attractiveParam);
	void spPosVerlet();
	void spVelVerlet(std::vector<double>& lenscales);
	void spNVE(std::vector<double>& lenscales, int nt);


	// Zebrafish active flow functions
	void dpInitializeActiveZebrafish(int NV, double phiDisk, double sizeDispersion, double R0);
	void dpActiveZebrafishNVE(double v0);
	void dpPntNVE(double v0);
	void dpZebrafishFireMinP(double Ptol, double Ktol);
	void dpPntFireMinP(double Ptol, double Ktol);
	void dpActiveZebrafishIsoCompression(double Ptol, double phiTarget, double deltaPhi);
	void dpPntIsoCompression(double Ptol, double phiTarget, double deltaPhi);
	void dpActiveZebrafishWallForces(double& wallPressure);
	void dpPNTWallForces();
	double zfishArea();
	double dpZfishPackingFraction();
	void dpZebrafishPositions();
	void dpPntPositions();
	void dpZebrafishEnergy();
	void dpPntEnergy();

	void initializeActiveZebrafish(std::vector<double>& radii, int NV, double phiDisk, double sizeDispersion, double R0);
	void fireMinimizeZebrafishSP(std::vector<double>& radii, double attractiveParam);
	void spActiveZebrafishWallForces(std::vector<double>& radii, double& wallPressure);
	void spActiveZebrafishNVE(std::vector<double>& radii, double v0);
	void spActiveZebrafishVicsek(std::vector<double>& radii, double attractionParam, double v0, double Dr, double vtau, double Pthresh, double dh);
	void printPositionsZebrafishSP(std::vector<double>& radii, std::vector<double>& psi);
	void printPositionsZebrafishSP(std::vector<double>& radii, std::vector<double>& psi, int Ncurr);
	void printEnergyZebrafishSP(double K);
	void printEnergyZebrafishSP(double K, int Ncurr);

	// active pipeflow functions
	void initializeActiveStickySP(std::vector<double>& radii, int NV, double phiDisk, double sizeDispersion, double Lscale);
	void pipeFireMinP(std::vector<double>& radii, double attractiveParam);
	void singleActiveCell(int NV, double phiInit, double calA0, double Dc, double vari, double v0);
	void spActiveForces(std::vector<double>& radii);
	void spActivePipeWallForces(std::vector<double>& radii);
	void spActivePipeNVE(std::vector<double>& radii, double T0);
	void spActivePipeFlow(std::vector<double>& radii, double a, double v0, double Dr);

	// Hopper functions
	void initializeHopperSP(std::vector<double>& radii, double w0, double w, double th, double Lmin);
	void initializeHopperDP(std::vector<double>& radii, double w0, double w, double th, double Lmin, int NV);
	void fireMinimizeHopperSP(std::vector<double>& radii, double w0, double w, double th);
	void fireMinimizeHopperDP(double w0, double w, double th);
	void fireMinimizeHopperF(double w0, double w, double th, double g);
	void hopperForcesSP(std::vector<double>& radii, double w0, double w, double th, double g, int closed);
	virtual void hopperForces(double w0, double w, double th, double g, int closed);
	void hopperWallForcesSP(std::vector<double>& radii, double w0, double w, double th, int closed);
	void hopperWallForcesDP(double w0, double w, double th, int closed);
	void hopperSPNVE(std::vector<double>& radii, double w0, double w, double th, double T0);
	void hopperDPNVE(double w0, double w, double th, double g, double T0);
	void flowHopperSP(std::vector<double>& radii, double w0, double w, double th, double g);
	int flowHopperDP(double w0, double w, double th, double g, double b);
	double hopperPackingFraction(std::vector<double>& radii, double w0, double w, double th);
	void setHopperPackingFraction(std::vector<double>& radii, double phiNew, double w0, double w, double th);
	void hopperPosVerletSP();
	void hopperVelVerletSP(std::vector<double>& radii, double b);
	void printHopperSP(std::vector<double>& radii, double w0, double w, double th, double g);
	void printHopperDP(double w0, double w, double th);
	double vertexWallForce(int onTheLeftC, int onTheLeftV, std::vector<double> & posX, std::vector<double> & posY);
	void hopperWallFriction(double w0, double w, double th, int closed);

	// printers
	void printSystemPositions();
	void printJammedConfig();
	void printSystemEnergy();
	void printSystemEnergy(int kmin);
	void printSystemContacts();
	void printSystemPositions(int frame);
	void printSystemEnergy(int frame, double Uval, double Kval);
	void printSystemStats();

	// activity

	//void activity(double T, double v0, double Dr, double vtau);
	void activityCOM(double T, double v0, double Dr, double vtau, double t_scale);
	void activityCOM_brownian(double T, double v0, double Dr, double vtau, double t_scale, int frames);
	void activityCOM_brownian_test(double T, double v0, double Dr, double vtau, double t_scale);
	void printCalA();
	void printContact();
	void printV();
	void conserve_momentum();
	void conserve_momentum_probe();
	double scale_v(double v0);

	void cell_NVE(double T, double v0, double Dr, double vtau, double t_scale, int frames);
	void random_shape_NVE(double T, double v0, double Dr, double vtau, double t_scale, int frames);
	void cell_NVE_probe(double T, double v0, double Dr, double vtau, double t_scale, int frames);
	void sp_NVE(double T, double v0, double Dr, double vtau, double t_scale, int frames);
	double * sp_NVE_tao(double T, double v0, double Dr, double vtau, double t_scale, int frames);
	double* NVE_tao(double T, double v0, double Dr, double vtau, double t_scale, int frames);
	int calTao(double q, int frames, std::vector< std::vector<double>>& x_com, std::vector< std::vector<double>>& y_com);
	void sp_NVE_probe(double T, double v0, double Dr, double vtau, double t_scale, int frames);
	void add_drag(int index, double force);
	void add_drag_cell(int index, double force);
	void sp_VelVerlet();
	void sp_VelVerlet(double b);
	void sp_VelVerlet_Langevin(double drag, double KbT, std::normal_distribution<double> & dist, std::mt19937 & gen);
	void sp_Forces(vector<double>& lenscales);
	void sp_Forces_probe(vector<double>& lenscales);
	void bumpy_NVE(double phiTarget, double T, double v0, double Dr, double vtau, double t_scale, int frames);
	virtual void bumpy_Forces();
	void bumpyRotation();
	void bumpy_angularV();
	void bumpy_angularV(double b);
	double totalRotaionalK();
	double cal_temp(double scaled_v);
	void rescal_V(double E);
	void rescal_V_probe(double E);




	void relaxP(double Ktolerance, double Ptolerance) {
		fireMinimizeP(Ptolerance, Ktolerance);
		phi = packingFraction();
		printJammedConfig_yc();
		phiPrintObject << phi << std::endl;
		printCalA();
		printContact();
	};
	void relaxF(double Ftolerance, double Ptolerance) {
		double Fckeck, Kcheck;
		fireMinimizeF(Ftolerance, Fckeck, Kcheck);
		phi = packingFraction();
		printJammedConfig_yc();
		phiPrintObject << phi << std::endl;
		printCalA();
		printContact();
	};
	void closeF() {
		jamPrintObject.close();
		packingPrintObject.close();
		energyPrintObject.close();
		statPrintObject.close();
		lengthscalePrintObject.close();
		phiPrintObject.close();
		calAPrintObject.close();
		contactPrintObject.close();
		vPrintObject.close();
	};
	void printJammedConfig_yc() {
		for (int ci = 0; ci < NCELLS; ci++)
			cell(ci).printVertexPositions_yc(jamPrintObject, ci);
	};
};

#endif