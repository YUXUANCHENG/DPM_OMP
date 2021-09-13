#ifndef CLI_HOPPER_H
#define CLI_HOPPER_H

#include "CLI.h"
#include <iomanip>
#include <iostream>
#include <vector>
#include <nlopt.hpp>

template <class Ptype = cellPacking2D>
class DPM_Hopper_CLI : public DPM_CLI<Ptype> {
public:

	double smallRadius = 0.5;			// radius fo smaller particles (diameter is length unit)
	double sizeRatio = 1.4;			// ratio of small diameter to large diameter
	//const double sizeRatio = 1;
	double w0 = 20.0;			// width of hopper reservoir (in units of small diameter)

	double th = PI / 4.0;		// hopper angle (pi - th = deflection angle from horizontal)
	double phi0 = 0.4;			// initial packing fraction
	double b = 0.1;
	//double b = 0.01;
	double g = 0.05;
	vector<double> radii;
	double w_scale = 2;
	double w;
	double gamafactor1 = 0;
	double gamafactor2 = 0;
	bool deformFlag = false;

	DPM_Hopper_CLI() {
		this->NT = 5e6;			// number of time steps for flow simulation
		this->NPRINT = 1e2;			// number of steps between printing
		// this->kl = 0.01;
		// this->ka = 0.5;
		this->kl = 1;
		this->ka = 10;
		this->kb = 0;
		// this->gamafactor1 = 0.9;
		// this->gamafactor2 = -1;
		this->g = 0.05;
		//this->Lini = 5;
		this->NBx = 30;
		this->NBy = 10;
		this->NCELLS = 800;
		this->NV = 16;
		this->calA0 = 1.0;
		this->kint = 2.0;
		this->Lini = this->NCELLS * (PI / 4) * (1 + sizeRatio * sizeRatio)/ 2/ 0.75 / pow(w0, 2);
		cout << "Lini = " << this->Lini << endl;
		//this->timeStepMag = 0.001;		
		this->radii = vector<double>(this->NCELLS, 0.0);
		for (int ci = 0; ci < this->NCELLS; ci++) {
			if (ci % 2 == 0)
				radii.at(ci) = smallRadius;
			else
				radii.at(ci) = smallRadius * sizeRatio;
		}	
	}

	virtual void setIndex(char const* argv[]) {
		string index_i_str = argv[1];
		string index_j_str = argv[2];
		stringstream index_i_ss(index_i_str);
		stringstream index_j_ss(index_j_str);
		index_i_ss >> this->index_i;
		index_j_ss >> this->index_j;
	}

	virtual void setKB() {
		;
	}

	virtual void setSeed() {
		//seed = index_i;
		this->seed = 1;
	}

	virtual void prepareSystem() {
		//w_scale = 0.5 + 0.05 * this->index_j;
		w_scale = 0.5 + 0.15 * this->index_j;
		w = w_scale * (1 + sizeRatio) / 2;
		// Initialze the system as disks
		cout << "	** Initializing hopper " << endl;
		if (deformFlag)
			this->particles->gDire = 1;
		this->particles->initializeHopperDP(radii, w0, w, th, this->Lini, this->NV);
		this->particles->forceVals(this->calA0, this->kl, this->ka, this->gam, this->kb, this->kint, this->del, this->aInitial);
		this->particles->initialize_subsystems(this->NBx, this->NBy);
		this->particles->vertexDPMTimeScale(this->timeStepMag);
		this->particles->changeL0(gamafactor1, gamafactor2);
		this->particles->closeF();
	}

	virtual void hopperFlow(char const* argv[])
	{
		this->qscompress(argv);
		_hopperFlow();
	}

	virtual void deformation(char const* argv[])
	{
		deformFlag = true;
		this->NT = 1e7;
		this->th = PI/4;
		//this->timeStepMag = 0.00002;	
		this->timeStepMag =  0.0001;	
		//this->timeStepMag = 0.0005;	
		this->kl = 1.74307;
		this->ka = 1;
		this->kb = 0;
		this->gamafactor1 = 0.0716771;
		this->gamafactor2 = -0.0345413;
		this->g = 0.170959;
		this->NCELLS = 1;
		this->NV = 64;
		this->w0 = 3;
		this->Lini = 5;
		this->qscompress(argv);
		this->particles->gDire = 1;
		this->particles->gOn = 0;
		_hopperFlow();
		double originalHeight = this->particles->calOriginalHeight();
		this->particles->gOn = 1;
		this->particles->hopperSimulation(w0, w, th, g, b);
		double endHeight = this->particles->calHeight();
		double angle = this->particles->calContactAng();
		double length = this->particles->calContactLength();
		cout << angle << endl;
		std::ofstream deformationPrint;
		deformationPrint.open("deformation.txt");
		deformationPrint << (originalHeight - endHeight)/originalHeight << "," << angle << "," << length <<endl;
		deformationPrint.close();

	}

	//void findParameter(char const* argv[]){
	double findParameter(const std::vector<double> &x){
		double presetCalA = 1.093;
		double presetAngle = 20;
		double threashold = 0.3;

		deformFlag = true;
		this->NT = 5e6;
		this->NPRINT = 1e3;	
		this->th = PI/4;
		this->g = x[3];
		this->timeStepMag = 0.0005;	
		this->kl = x[0];
		this->gamafactor1 = x[1];
		this->gamafactor2 = x[2];
		this->ka = 1;
		this->kb = 0;
		this->NCELLS = 1;
		this->NV = 64;
		this->w0 = 3;
		this->Lini = 5;
		const char *right[3] = {"2", "0", "10"};
		this->qscompress(right);
		this->particles->gDire = 1;
		this->extend = "_" + to_string(this->index_i) + to_string(this->index_j) + ".txt";
		string energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF;
		this->produceFileName(this->extend, energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF);
		this->particles->openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF);
		this->particles->gOn = 1;
		this->particles->fireMinimizeHopperF(w0, w, th, g);
		std::ofstream debugInfo;
		debugInfo.open("debug.txt", std::ios_base::app);
		debugInfo << "kl = " << this->kl << ", gama1 = " << this->gamafactor1 << ", gama2 = " << this->gamafactor2 << ", g = " << this->g << endl;
		return this->particles->matchPreset(presetCalA, presetAngle, threashold);
					
	}

	virtual void _hopperFlow() {
		this->extend = "_" + to_string(this->index_i) + to_string(this->index_j) + ".txt";
		string energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF;
		this->produceFileName(this->extend, energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF);
		this->particles->openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF);
		//this->particles->initialize_subsystems();
		int result = this->particles->hopperSimulation(w0, w, th, g, b);
		this->v0PrintObject << this->kl << "," << this->gam << "," << g << "," << w_scale << "," << result << "," << this->kb << endl;
		cout << "	** FINISHED **   " << endl;
		//clogPrintObject.close();
		this->v0PrintObject.close();
	}
};

template <class Ptype = Bumpy>
class Bumpy_Hopper_CLI : public DPM_Hopper_CLI<Ptype> {
public:
	//typedef Bumpy particleType;
	Bumpy_Hopper_CLI() : DPM_Hopper_CLI<Ptype>() {
	this->kl = 0;
	this->ka = 0;
	this->kb = 0;
	}
	virtual void prepareSystem() {
		DPM_Hopper_CLI<Ptype>::prepareSystem();
		this->particles->calInertia();
	}
/*
	virtual void createParticles(char const* argv[])
	{
		_createParticles<Bumpy>(argv);
	}
*/
};

double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    if (!grad.empty()) {
        // grad[0] = 0.0;
        // grad[1] = 0.5 / sqrt(x[1]);
		;
    }
	DPM_Hopper_CLI<> cli;
	return cli.findParameter(x);
    //return sqrt(x[1]);
}

// typedef struct {
//     double kl, gama1, gama2;
// } my_constraint_data;

// double myconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data)
// {
//     my_constraint_data *d = reinterpret_cast<my_constraint_data*>(data);
//     double a = d->a, b = d->b;
//     if (!grad.empty()) {
//         grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
//         grad[1] = -1.0;
//     }
//     return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
// }

void optFun()
{
	nlopt::opt opt(nlopt::LN_COBYLA, 4);
	std::vector<double> lb(4);
	lb[0] = 0.1; lb[1] = 0.01; lb[2] = -0.1; lb[3] = 0.05;
	opt.set_lower_bounds(lb);
	std::vector<double> ub(4);
	ub[0] = 3; ub[1] = 0.2; ub[2] = 0.5; ub[3] = 0.6;
	opt.set_upper_bounds(ub);
	opt.set_min_objective(myfunc, NULL);
	// my_constraint_data data[2] = { {2,0}, {-1,1} };
	// opt.add_inequality_constraint(myconstraint, &data[0], 1e-8);
	// opt.add_inequality_constraint(myconstraint, &data[1], 1e-8);
	opt.set_xtol_rel(1e-4);
	std::vector<double> x(4);
	x[0] = 1; x[1] = 0.1; x[2] = 0.3; x[3] = 0.3;
	double minf;

	try{
		nlopt::result result = opt.optimize(x, minf);
		std::cout << "found minimum at f(" << x[0] << "," << x[1] << "," << x[2] << "," << x[3] <<") = "
			<< std::setprecision(5) << minf << std::endl;
	}
	catch(std::exception &e) {
		std::cout << "nlopt failed: " << e.what() << std::endl;
	}
}

#endif