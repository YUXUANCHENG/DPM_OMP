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

#include "deformableParticles2D.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

class cellPacking2D{
private:

	// int scalars
	int NDIM;						// spatial dimension (will always be 2)
	int NCELLS;						// number of cells
	int NT;							// number of total time steps
	int NPRINT;						// number of time steps between print steps

	// double scalars
	double seed;					// initial seed
	double dt;						// time step size
	double dt0;						// initial time step size (for dynamic time step methods)
	double T;						// system temperature
	double phi;						// system packing/volume fraction

	// boundary lengths
	std::vector<double> L;

	// virial stresses
	double sigmaXX, sigmaXY, sigmaYX, sigmaYY;
	
	// array of cells
	deformableParticles2D* cellArray;

	// contact matrix
	int* contactMatrix;				// array of contacts between particles

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

public:

	// Constructors and Destructors
	void defaultvars();
	cellPacking2D();
	cellPacking2D(int ncells, int nt, int nprint, double l, double s);
	cellPacking2D(int ncells, int ntumor, int tumorNV, int adiposeNV, double tumorCalA, double adiposeCalA, int s);
	cellPacking2D(std::ifstream& inputFileObject, double asphericity, double s);
	~cellPacking2D();

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
	
	void openJamObject(std::string& str, std::string& str1, std::string& str2, std::string& str3, std::string& str4, std::string& str5) {
		jamPrintObject.open(str.c_str());
		lengthscalePrintObject.open(str1.c_str());
		phiPrintObject.open(str2.c_str());
		calAPrintObject.open(str3.c_str());
		contactPrintObject.open(str4.c_str());
		vPrintObject.open(str5.c_str());
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

	// box len
	double getL(int d) { return L.at(d); };

	deformableParticles2D& cell(int ci);	// return cell object ci
	int nframes();							// number of frames in the simulation
	int cmindex(int ci, int cj);			// contact matrix index
	int contacts(int ci, int cj);			// contact matrix element
	int totalContacts();

	// system-wide calculations
	int totalNumberOfContacts();			// calculate total number of contacts
	double packingFraction();				// calculate system packing fraction
	double totalPotentialEnergy();			// calculate total potential energy
	double interactionPotentialEnergy();	// calculate only interaction potential energy
	double totalKineticEnergy();			// calculate total kinetic energy
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

	// set force values for all cells to be the same
	void forceVals(double calA0, double kl, double ka, double gam, double kb, double kint, double del, double a);

	
	void vertexDPMTimeScale(double timeStepMag);
	void addContact(int ci, int cj);
	void deleteContact(int ci, int cj);
	void resetContacts();
	int particleContacts(int ci);
	void setPackingFraction(double val);
	void setAsphericity(double val);
	void setAsphericity(int ci, double val);
	void scaleLengths(double val);
	void rescaleVelocities(double temperature);
	int removeRattlers(int krcrs);

	/**************************

		Forces and position 
			updates

	***************************/

	void calculateForces();
	void gelationForces();


	/**************************

		FIRE energy minimzation

	***************************/

	// FIRE 2.0 relaxation functions
	void fireMinimizeP(double Ptol, double Ktol);
	void fireMinimizeF(double Ftol, double Ktol, double& Ftest, double& Ktest);


	/**************************

		Vibrational 
			density of states

	***************************/

	void cellVDOS(std::ofstream& vdosOutObj, double dphi, double Ftol, double Ktol);


	/**************************

		Simulation Functions:

	***************************/

	// NVE test function
	void cellNVE();
	void cellOverDamped();

	// Find jammed state from initially dilute configuration
	void findJamming(double dphi0, double Ktol, double Ftol, double Ptol);

	// compress isotropically to a target packing fraction
	void qsIsoCompression(double phiTarget, double deltaPhi, double Ftol, double Ktol);

	// Gelation functions
	void twoParticleContact(int NV);
	void initializeGel(int NV, double phiDisk, double sizeDispersion, double delval);
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
	void initializeHopperSP(std::vector<double>& radii, double w0, double w, double th, double Lmin, int NV);
	void fireMinimizeHopperSP(std::vector<double>& radii, double w0, double w, double th);
	void fireMinimizeHopperDP(double w0, double w, double th);
	void hopperForcesSP(std::vector<double>& radii, double w0, double w, double th, double g, int closed);
	void hopperForcesDP(double w0, double w, double th, double g, int closed);
	void hopperWallForcesSP(std::vector<double>& radii, double w0, double w, double th, int closed);
	void hopperWallForcesDP(double w0, double w, double th, int closed);
	void hopperSPNVE(std::vector<double>& radii, double w0, double w, double th, double T0);
	void flowHopperSP(std::vector<double>& radii, double w0, double w, double th, double g);
	void flowHopperDP(double w0, double w, double th, double g);
	double hopperPackingFraction(std::vector<double>& radii, double w0, double w, double th);
	void setHopperPackingFraction(std::vector<double>& radii, double phiNew, double w0, double w, double th);
	void hopperPosVerletSP();
	void hopperVelVerletSP(std::vector<double>& radii);
	void printHopperSP(std::vector<double>& radii, double w0, double w, double th, double g);
	void printHopperDP(double w0, double w, double th, double g);

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
	double scale_v(double v0);


	void relaxP(double Ktolerance, double Ptolerance) {
		fireMinimizeP(Ptolerance, Ktolerance);
		phi = packingFraction();
		printJammedConfig_yc();
		phiPrintObject << phi << std::endl;
		printCalA();
		printContact();
	};
	void relaxF(double Ktolerance, double Ftolerance, double Ptolerance) {
		double Fckeck, Kcheck;
		fireMinimizeF(Ftolerance, Ktolerance, Fckeck, Kcheck);
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