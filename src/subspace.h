#ifndef SUBSPACE_H
#define SUBSPACE_H

#include "deformableParticles2D.h"
#include "cellPacking2D.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <stack>

class cvpair {
public:
	int ci;
	int vi;
	int boxid = 0;

	cvpair(int c, int v) {
		ci = c;
		vi = v;
	}
};
using namespace std;

class cellPacking2D;

class subspace {

protected:

	// resident cells and cashed cells
	vector<cvpair*> resident_cells;
	vector<cvpair*> cashed_cells;

	// list indicates near boundary cells that need to be sent to neighbor boxes
	vector<cvpair*> cash_out_list_lower;
	vector<cvpair*> cash_out_list_upper;

	// pointer to the whole system (cell_group)
	cellPacking2D* pointer_to_system;

	// which box is this
	int box_id;

	vector<double> L;
	vector<int> N_systems;

	int NDIM = 2;
	double Ncc = 0, Nvv = 0;

	// stress
	double sigmaXX = 0.0;
	double sigmaXY = 0.0;
	double sigmaYX = 0.0;
	double sigmaYY = 0.0;

	// MD time 
	double dt0;
	double PI = 4 * atan(1);

	// Update frequency
	// seems like cashe frequency won't affect speed too much
	int update_freqency = 10;

	// indicate what fraction of the system size will be cashed
	vector<double> cashed_fraction{ 0.1, 0.1 };
	double cashed_length = 2;

public:

	void initialize(cellPacking2D* const& pointer, vector<double> const& L, vector<int> const& N_systems, int box_id, double const& dt0) {
		pointer_to_system = pointer;
		this->L = L;
		this->box_id = box_id;
		this->N_systems = N_systems;
		this->dt0 = dt0;
	};

	void cashe_out(int direction);
	void reset_cashe();
	void reset();
	int neighbor_box(int direction, int upper_lower);
	double find_boundary(int direction, int upper_lower);
	void migrate_out();

	void cashe_in(vector<cvpair*>& cash_list);
	void migrate_in(cvpair* const& migration);

	void calculateForces_insub();
	void calculateForces_betweensub();
	void activityCOM_brownian_insub(double T, double v0, double Dr, double vtau, double t_scale, int frames);
	void fireMinimizeF_insub(double Ftol, double& Fcheck, double& Kcheck, double& P, double& vstarnrm, double& fstarnrm, bool& converged);
	double forceRMS_insub();
	double totalKineticEnergy_insub();
	double max_length();
	void print_information();
	void cal_cashed_fraction();
	int vertexForce(cvpair* onTheLeft, cvpair* onTheRight, double& sigmaXX, double& sigmaXY, double& sigmaYX, double& sigmaYY);
	

};







#endif