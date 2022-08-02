#ifndef SUBSPACE_H
#define SUBSPACE_H

#include "deformableParticles2D.h"
//#include "DPM_Parallel.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <stack>
#include <Eigen/Core>
#include <Eigen/Dense>
extern bool settleDown;

typedef Eigen::Matrix<double, 2,  2>  MATRIX2;
typedef Eigen::Matrix<double, 3,  3>  MATRIX3;
typedef Eigen::Matrix<double, 4,  4>  MATRIX4;
typedef Eigen::Matrix<double, 6,  6>  MATRIX6;
typedef Eigen::Matrix<double, 9,  9>  MATRIX9;
typedef Eigen::Matrix<double, 4,  6>  MATRIX4x6;
typedef Eigen::Matrix<double, 3,  12> MATRIX3x12;
typedef Eigen::Matrix<double, 9,  12> MATRIX9x12;
typedef Eigen::Matrix<double, 12, 12> MATRIX12;
typedef Eigen::Matrix<double, 2,  1>  VECTOR2;
typedef Eigen::Matrix<double, 3,  1>  VECTOR3;
typedef Eigen::Matrix<double, 4,  1>  VECTOR4;
typedef Eigen::Matrix<double, 6,  1>  VECTOR6;
typedef Eigen::Matrix<double, 9,  1>  VECTOR9;
typedef Eigen::Matrix<double, 12, 1>  VECTOR12;
typedef Eigen::Matrix<double, 2,  6> MATRIX2x6;

typedef Eigen::Matrix<int, 2, 1> VECTOR2I;
typedef Eigen::Matrix<int, 3, 1> VECTOR3I;
typedef Eigen::Matrix<int, 4, 1> VECTOR4I;

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VECTOR;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MATRIX;

extern bool frictionFlag;


class cvpair {
public:
	int ci;
	int vi;
	int boxid = 0;

	cvpair(int c, int v) {
		ci = c;
		vi = v;
	}
	cvpair(){ci = -1; vi = -1;};
	bool operator<(const cvpair& src)const
    {
        return (this->ci < src.ci || this->vi < src.vi);
    }
};
using namespace std;

class DPM_Parallel;

class subspace {
	
public:
	double Ncc = 0, Nvv = 0;

	// stress
	double sigmaXX = 0.0;
	double sigmaXY = 0.0;
	double sigmaYX = 0.0;
	double sigmaYY = 0.0;
	vector<double> L;
	double coefu = 1, coefV = 10;

	// resident cells and cashed cells
	vector<cvpair*> resident_cells;
	vector<cvpair*> cashed_cells;

	// list indicates near boundary cells that need to be sent to neighbor boxes
	vector<cvpair*> cash_out_list_lower;
	vector<cvpair*> cash_out_list_upper;

	// pointer to the whole system (cell_group)
	DPM_Parallel* pointer_to_system = nullptr;

	// which box is this
	int box_id;

	
	vector<int> N_systems;

	int NDIM = 2;

	// MD time 
	double dt0;

	// Update frequency
	// seems like cashe frequency won't affect speed too much
	int update_freqency = 10;

	// indicate what fraction of the system size will be cashed
	vector<double> cashed_fraction{ 0.1, 0.1 };
	double cashed_length = 2;


	void initialize(DPM_Parallel* const& pointer, vector<double> const& L, vector<int> const& N_systems, int box_id, double const& dt0) {
		pointer_to_system = pointer;
		this->L.reserve(2);
		this->L = L;
		this->box_id = box_id;
		this->N_systems = N_systems;
		this->dt0 = dt0;
	};

	void setFriction(double coefu, double coefv){
		this->coefu = coefu;
		this->coefV = coefv;
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
	double totalKineticEnergy_insub();
	void print_information();
	void cal_cashed_fraction();
	int vertexForce(cvpair* onTheLeft, cvpair* onTheRight, double& sigmaXX, double& sigmaXY, double& sigmaYX, double& sigmaYY);
	int vertexForce_with_Torque(cvpair* onTheLeft, cvpair* onTheRight, double& sigmaXX, double& sigmaXY, double& sigmaYX, double& sigmaYY);
	virtual double vertexEdgeDist(const cvpair* onTheLeft, const cvpair* onTheRight){return 0;};

};

class frictionlessSubspace: public virtual subspace {
public:

	void calculateEdgeForces_insub();
	void calculateEdgeForces_betweensub();
	int vertexEdgeForce(cvpair* onTheLeft, cvpair* onTheRight, double& sigmaXX, double& sigmaXY, double& sigmaYX, double& sigmaYY);
	virtual double vertexEdgeDist(const cvpair* onTheLeft, const cvpair* onTheRight);
};





#endif