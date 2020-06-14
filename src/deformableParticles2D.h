#ifndef DEFORMABLEPARTICLES2D_H
#define DEFORMABLEPARTICLES2D_H


/*

	-- DEFORMABLE PARTICLE CLASS -- 

	Class for individual deformable particles, 
	for use as member variables in an MD class

	Holds position, velocity and force information

	note: everything in d=2, not generalizable to any d

*/

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
using namespace std;

class deformableParticles2D{
private:
	// spatial dimension
	int NDIM;

	// number of vertices
	int NV;

	// periodic boundary conditions
	std::vector<int> pbc;

	// box lengths (Lx and Ly)
	std::vector<double> L;

	// vertex positions relative to center of mass
	double* vertexPositions;

	// vertex velocity in lab frame
	double* vertexVelocity;

	// vertex accerlation in lab frame
	double* vertexAcceleration;

	// forces on vertices
	double* vertexForces;

	// cell position in lab frame
	double* cellPosition;

	// interaction potential energy
	double* interactionPotential;

	// force parameters
	double kl;			// perimeter energy scale
	double ka;			// area energy scale
	double gam;			// surface tension energy scale
	double kb;			// bending energy scale
	double kint;		// interaction energy scale

	// rest parameters
	double l0;			// rest length for each perimeter spring
	double a0;			// rest area for particles
	double c0; 			// rest angle for cell boundary
	double del;			// contact distance for two edge segments
	double a;			// only attraction parameter (distance and strength, good for DM)

	double c_psi = 2 * 3.1415 * (double)rand() / (RAND_MAX + 1.0);

	// shear strain (for LEbc)
	double strain;

public:
	// constructors
	deformableParticles2D();
	deformableParticles2D(int n);
	~deformableParticles2D();

	// operators
	void operator=(deformableParticles2D& onTheRight);

	// initialization
	void initializeVertices();
	void initializeCell();
	void regularPolygon();
	void vertexPerturbation(double dscale);

	// getters (simple access)
	int getNV() { return NV; };
	double getkl() { return kl; };
	double getka() { return ka; };
	double getgam() { return gam; };
	double getkb() { return kb; };
	double getkint() { return kint; };
	double getl0() { return l0; };
	double geta0() { return a0; };
	double getc0() { return c0; };
	double getdel() { return del; };
	double geta() { return a; };
	double getstrain() { return strain; };

	// access pbc and box length information
	int getpbc(int d) { return pbc.at(d); };
	double getL(int d) { return L.at(d); };

	// getters (defined in .cpp file)
	double vpos(int vertex, int dim);
	double vrel(int vertex, int dim);
	double vvel(int vertex, int dim);
	double vacc(int vertex, int dim);
	double vforce(int vertex, int dim);
	double cpos(int dim);
	double cvel(int dim);
	double cforce(int dim);
	double uInt(int vertex);
	double distance(deformableParticles2D& onTheRight, int vj, int vi, int d);
	double cellDistance(deformableParticles2D& onTheRight, int d);

	// setters (simple mutation)
	void setNV(int nv) { NV = nv; };
	void setkl(double val) { kl = val; };
	void setka(double val) { ka = val; };
	void setgam(double val) { gam = val; };
	void setkb(double val) { kb = val; };
	void setkint(double val) { kint = val; };
	void setl0(double val) { l0 = val; };
	void seta0(double val) { a0 = val; };
	void setc0(double val) { c0 = val; };
	void setc0Angle(double val) { c0 = cos(val); };
	void setdel(double val) { del = val; };
	void seta(double val) { a = val; };
	void setpbc(int d, int val) { pbc.at(d) = val; };
	void setL(int d, double val) { L.at(d) = val; };
	void setstrain(double val) { strain = val; };

	// setters (defined in .cpp file)
	void setVPos(int vertex, int dim, double val);
	void setVRel(int vertex, int dim, double val);
	void setVVel(int vertex, int dim, double val);
	void setVAcc(int vertex, int dim, double val);
	void setVForce(int vertex, int dim, double val);
	void setCPos(int dim, double val);
	void setCVel(int dim, double val);
	void setCForce(int dim, double val);
	void setUInt(int vertex, double val);
	void setAsphericity(double val);
	void setAsphericityConstA(double val);

	// update cpos
	void updateCPos();

	// scale all lengths
	void scale(double val);

	// calculations
	double polygonArea();							// area of polygon (without vertices)
	double safe_acos(double x);
	double area();									// area of cell
	double perimeter();								// perimeter of cell
	double asphericity();							// instantaneous asphericity parameter
	double calA0();									// preferred asphericity
	double segmentLength(int vertex);				// distance between vertex i and i+1
	double segment(int vertex, int dim);			// vector connecting i and i+1
	double dotProduct(int v1, int v2);				// dot product between vertices v1 and v2
	double segmentDotProduct(int l1, int l2);		// dot product between segments l1 and l2
	double segmentCosine(int l1);					// dot product between segment l1 - 1 and l1

	// force functions
	void shapeForces();
	int segmentForce(deformableParticles2D& onTheRight); // return 0 or 1, depending on contact
	int vertexForce(deformableParticles2D& onTheRight, double& sigmaXX, double& sigmaXY, double& sigmaYX, double& sigmaYY);
	int vertexForce(deformableParticles2D &onTheRight, double& sigmaXX, double& sigmaXY, double& sigmaYX, double& sigmaYY, double aij);
	int pwAttractiveContacts(deformableParticles2D &onTheRight);
	int radialForce(deformableParticles2D& onTheRight, double bscale); 

	// energy functions
	double perimeterEnergy();
	double areaEnergy();
	double surfaceTensionEnergy();
	double bendEnergy();
	double interactionEnergy();
	double totalPotentialEnergy();
	double totalKineticEnergy();

	// overloaded operators
	void operator% (const deformableParticles2D &onTheRight); 	// interaction force calculation operator

	// integrator options
	void verletPositionUpdate(double dt);
	void BrownianPositionUpdate(double dt);
	void verletVelocityUpdate(double dt);
	void verletVelocityUpdate(double dt, double dampingParam);
	
	// print functions
	void printVertexPositions(std::ofstream& vertexPrintObject, int cellID);
	void printVertexPositions(std::ofstream& vertexPrintObject, int cellID, int frame);
	void printCellEnergy(std::ofstream& energyPrintObject, int frame);

	void activeVerletVelocityUpdateCOM(double dt0, double Dr, double vtau, double v0);
	void activeVerletVelocityUpdateCOM_brownian(double dt0, double Dr, double vtau, double v0);

	void printVertexPositions_yc(std::ofstream& vertexPrintObject, int cellID) {
		int p = 14;
		// loop over vertices, print
		for (int i = 0; i < NV; i++) {
			for (int d = 0; d < NDIM; d++) {
				/*
				if (vpos(i, d) - cpos(d) > 0.5 * L.at(d))
					vertexPrintObject << std::setprecision(p) << vpos(i, d) - L.at(d);
				else if (cpos(d) - vpos(i, d) > 0.5 * L.at(d))
					vertexPrintObject << std::setprecision(p) << vpos(i, d) + L.at(d);
				else
				*/
					vertexPrintObject << std::setprecision(p) << vpos(i, d);

				vertexPrintObject << ",";
			}
			vertexPrintObject << l0 << std::endl;


		}
	};

	double calA();
	void printlengthscale(std::ofstream& lengthscalePrintObject) {

		lengthscalePrintObject << NV << std::endl;

	};

	double cal_mean_v(int d);

	double momentum(int d);
};




#endif