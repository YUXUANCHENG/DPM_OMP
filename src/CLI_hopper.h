#ifndef CLI_HOPPER_H
#define CLI_HOPPER_H

#include "CLI.h"
#include <iomanip>
#include <iostream>
#include <vector>

#define optimizer 0

# if optimizer
#include <nlopt.hpp>
#endif
extern bool constPressureFlag;
extern bool frictionFlag;
extern bool replaceFlag;

// check list: kb ksp interval angle, mu factor, script, softflag, cpu, partition
template <class Ptype = cellPacking2D>
class DPM_Hopper_CLI : public DPM_CLI<Ptype> {
public:

	double smallRadius = 0.5;			// radius fo smaller particles (diameter is length unit)
	//double smallRadius = 0.3;			// radius fo smaller particles (diameter is length unit)
	double sizeRatio = 1.4;			// ratio of small diameter to large diameter
	//double sizeRatio = 4;			// ratio of small diameter to large diameter
	// const double sizeRatio = 1.01;
	double w0 = 20.0;			// width of hopper reservoir (in units of small diameter)
	// double w0 = 60.0;			// width of hopper reservoir (in units of small diameter)

	double th = PI / 4.0;		// hopper angle (pi - th = deflection angle from horizontal)
	// double th =  (90.0 - 89.0)/180 * PI;		// hopper angle (pi - th = deflection angle from horizontal)
	double phi0 = 0.4;			// initial packing fraction
	double b = 0.1;
	// double b = 0;
	double g = 0.05;
	vector<double> radii;
	double w_scale = 2;
	double w;
	double gamafactor1 = 0;
	double gamafactor2 = 0;
	bool deformFlag = false;
	double coefu = 1;
	double coefv = 0;
	double spK = 10;
	double start = 1;
	double interval = 0.1;

	double scaleFactor = 1;

	DPM_Hopper_CLI() {
		this->NT = 5e6;			// number of time steps for flow simulation
		// this->NT = 1e6;			// number of time steps for flow simulation
		this->NPRINT = 1e4;			// number of steps between printing
		this->kl = 0.01*scaleFactor;
		this->ka = 100*scaleFactor;
		this->kb = 0.001*scaleFactor;
		if (frictionFlag)
			this->b = 0;
		this->b *= scaleFactor;
		int factorNx = 5;
		// this->gamafactor1 = 0.9;
		// this->gamafactor2 = -1;
		//this->g = 0.05;
		if (constPressureFlag)
		{
			this->g = 0.5*scaleFactor;
		}
		else
			// this->g = 0.05*scaleFactor;
			this->g = 0.3*scaleFactor;
		// this->NBx = 30;
		// this->NBy = 10;
		int mixed = 1;
		if (mixed)
		{
			this->b = 1;
			this->coefu = 0;
			this->coefv = 0;
			// this->coefv = 10;
			this->spK = 1000;
			// this->g = 0.01;
		}
		this->kint = 2.0*scaleFactor;
		if (replaceFlag)
		{
			// this->NT = 1e8;	
			// this->NPRINT = 1e7;	
			this->NT = 1e6;	
			this->NPRINT = 1e3;	
			// this->NPRINT = 1e5;	
			// this->NCELLS = 1600;
			// this->NCELLS = 400;
			this->NCELLS = 1;
			// this->timeStepMag = 0.002;
			// w0 = 60.0;
			w0 = 25.0;
			th =  (90.0 - 89.9)/180 * PI;
			// th = PI/4.0;
			// th = PI/3.0;
			// th = PI/6.0;
			// th =  (20.0)/180 * PI;
			// th =  (15.0)/180 * PI;
			this->kint = 10*scaleFactor;
			factorNx = 1;
		}
		else
		{
			// this->NCELLS = 64;
			// this->NPRINT = 1e2;
			// this->NCELLS = 800;
			// this->NCELLS = 100;
			this->NCELLS = 200;
			if (this->NCELLS > 100)
			{
				this->NT = 1e8;			// number of time steps for flow simulation
				this->NPRINT = 1e5;	
				// this->timeStepMag = 0.002;
				w0 = 40.0;
			}
		}
		
		// this->NV = 64;
		this->NV = 16;
		this->calA0 = 1.0;
		// this->calA0 = 1.15;
		this->NBy = 4 * round(w0/10) * this->NV/16;
		this->NBx = 2 * factorNx * ceil((this->NCELLS/64.0)) * pow(this->NV/16, 1) * ceil((this->NBy/30.0));
		// this->NBx = 5 * (this->NCELLS/64) * pow(this->NV/16, 2) / (this->NBy/30);
		if (this->NCELLS > 300)
			this->Lini = this->NCELLS * (PI / 4) * (1 + sizeRatio * sizeRatio)/ 2/ 0.6 / pow(w0, 2);
		else
			this->Lini = 1;
		cout << "Lini = " << this->Lini << endl;
		// if (this->kb > 9 || this->NV > 16)
		// if (this->kb > 9)
		// 	// this->timeStepMag = 0.001;		
		// 	this->timeStepMag = 0.002;		
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
		// ;
		// this->kl = 0.01 * scaleFactor * pow(this->index_i%10 + 1, 3);
		// this->kb = 0.001 * scaleFactor * pow(this->index_i/10 + 1, 3);

		// this->kl = 0.1 * scaleFactor * pow(10, (this->index_i%10)*3.0/9);
		// this->kb = 0.001 * scaleFactor * pow(10, (this->index_i/10)*4.0/9);
		// this->spK = 1 * scaleFactor * pow(10, (this->index_i%10)*3.0/9);
		// this->coefv = 0.1 * (this->index_i + 1);
		// setFriction();
		// uglyNdepend();
		// if (this->kb > 9)
		// 	// this->timeStepMag = 0.001;		
		// 	this->timeStepMag = 0.002;
	}

	void setFriction()
	{
		int i = this->index_i;
		vector<double> bVec{1e-6,1e-5,1e-4,1e-3,3e-3,5e-3,8e-3,2e-2,4e-2,6e-2,8e-2,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0};
		vector<double> vVec{10,10,10,10,10,10,10,10,10,10,10,10,5,1,0.1,1e-2,1e-3,0,10};
		this->b = bVec.at(i);
		this->coefv = vVec.at(i);
		this->coefu = 0;
		// this->coefv *= 5;
		// this->coefv /= 10;
		// this->coefv /= 100;
		// this->coefu = this->coefv/10;
	}

	void uglyNdepend()
	{
		this->NCELLS = 800 + this->index_i * 320;
		this->NBx = 1 * ceil((this->NCELLS/64.0)) * pow(this->NV/16, 1) * ceil((this->NBy/30.0));
		this->Lini = this->NCELLS * (PI / 4) * (1 + sizeRatio * sizeRatio)/ 2/ 0.6 / pow(w0, 2);
		cout << "Lini = " << this->Lini << endl;
		this->radii = vector<double>(this->NCELLS, 0.0);
		for (int ci = 0; ci < this->NCELLS; ci++) {
			if (ci % 2 == 0)
				radii.at(ci) = smallRadius;
			else
				radii.at(ci) = smallRadius * sizeRatio;
		}	
	}

	virtual void setSeed() {
		this->seed = this->index_i;
		// this->seed = 0;
	}

	virtual void prepareSystem() {
		// start = 1; interval = 0.6;
		start = 3; interval = 0.4;
		// start = 5; interval = 0.5;
		// start = 0.5; interval = 0.1;
		// w_scale = 0.5 + 0.05 * this->index_j;
		// w_scale = 3 + 0.1 * this->index_j;
		// w_scale = 0.5 + 0.1 * this->index_j;
		// w_scale = 3 + 1 * this->index_j;
		// w_scale = start + interval * this->index_j;
		// w_scale = 2.2;
		// w_scale = 0.9;
		w_scale = 10;
		// w_scale = 2.5;
		// w_scale = 0.5 + 1 * this->index_j;
		w = w_scale * (1 + sizeRatio) / 2;
		double ObX;
		// int cutoff = 15; double spacing = 0.2;
		// if (this->index_j <= cutoff)
		// 	ObX = (1 + spacing * this->index_j) * (1 + sizeRatio) / 2;
		// else
		// 	ObX = (1 + spacing * cutoff + 1 * (this->index_j-cutoff) ) * (1 + sizeRatio) / 2;
		// double ObR = 6 * (1 + sizeRatio) / 2;
		double ObR = 10;
		ObX = 0;
		// Initialze the system as disks
		cout << "	** Initializing hopper " << endl;
		if (deformFlag)
			this->particles->gDire = 1;
		int fricTmp = frictionFlag;
		frictionFlag = 0;
		this->particles->setObstacle(ObX,0,ObR);
		this->particles->initializeHopperDP(radii, w0, w, th, this->Lini, this->NV);
		frictionFlag = fricTmp;
		this->particles->forceVals(this->calA0, this->kl, this->ka, this->gam, this->kb, this->kint, this->del, this->aInitial);
		this->particles->setFriction(this->coefu, this->coefv);
		this->particles->setSPk(this->spK);
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
		this->convert2double(argv[3], this->kl);
		this->convert2double(argv[4], this->gamafactor1);
		this->convert2double(argv[5], this->gamafactor2);
		this->convert2double(argv[6], this->g);

		this->NT = 1e6;
		this->NPRINT = 1e3;	
		// this->timeStepMag = 0.00005;	
		this->timeStepMag =  0.0005;	
		this->ka = 10*scaleFactor;
		this->kb = 10*scaleFactor;
		this->kl = this->kl*scaleFactor;
		this->g = this->g*scaleFactor;
		this->NCELLS = 1;
		this->NV = 16;
		this->w0 = 10;
		this->Lini = 5;
		// this->calA0 = 1.15;
		this->calA0 = 1;
		this->kint = 2.0*scaleFactor;
		this->qscompress(argv);

		this->particles->gDire = 1;
		this->extend = "_" + to_string(this->index_i) + to_string(this->index_j) + ".txt";
		string energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF;
		this->produceFileName(this->extend, energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF);
		this->particles->openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF);
		double originalHeight;
		// if (this->calA0 < 1 + 1e-7)
		if (false)
		{
			this->particles->gOn = 0;
			this->particles->fireMinimizeHopperF(w0, w, th, g);
			originalHeight = this->particles->calOriginalHeight();
		}
		else
		{
			this->particles->gOn = 1;
			this->particles->fireMinimizeHopperF(w0, w, th, 1e-5);
			originalHeight = this->particles->calHeight();
		}

		this->particles->gOn = 1;
		this->particles->fireMinimizeHopperF(w0, w, th, g);
		double endHeight = this->particles->calHeight();
		double angle = this->particles->calContactAng();
		double length = this->particles->calContactLength();
		cout << angle << endl;
		std::ofstream deformationPrint;
		deformationPrint.open("deformation.txt");
		deformationPrint << (originalHeight - endHeight)/originalHeight << "," << angle << "," << length <<endl;
		deformationPrint.close();

	}

	virtual void measureFriction(char const* argv[])
	{
		deformFlag = true;
		this->convert2double(argv[3], this->kl);
		this->convert2double(argv[4], this->gamafactor1);
		this->convert2double(argv[5], this->gamafactor2);
		this->convert2double(argv[6], this->g);

		this->NT = 1e6;
		this->NPRINT = 1e2;	
		// this->timeStepMag = 0.00005;	
		this->timeStepMag =  0.0005;	
		this->ka = 10*scaleFactor;
		this->kb = 0*scaleFactor;
		this->kl = this->kl*scaleFactor;
		this->g = this->g*scaleFactor;
		this->NCELLS = 1;
		this->NV = 64;
		this->w0 = 10;
		this->Lini = 5;
		// this->calA0 = 1.15;
		this->calA0 = 1;
		this->kint = 2.0*scaleFactor;
		this->qscompress(argv);

		this->particles->gDire = 1;
		this->particles->gOn = 0;

		this->particles->fireMinimizeHopperF(w0, w, th, g);


		this->particles->gDire = 0;
		this->particles->gOn = 1;
		_hopperFlow();
	}

	//void findParameter(char const* argv[]){
	double findParameter(const std::vector<double> &x){
		//double presetCalA = 1.142;
		//double presetAngle = 17.9;
		double presetCalA = 1.131;
		double presetAngle = 22.13;
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
		if (constPressureFlag)
			this->particles->gOn = 0;
		this->extend = "_" + to_string(this->index_i) + to_string(this->index_j) + ".txt";
		string energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF;
		this->produceFileName(this->extend, energyF, jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF);
		this->particles->openJamObject(jammingF, lengthscaleF, phiF, calAF, contactF, vF, ISF);
		//this->particles->initialize_subsystems();
		int result = this->particles->hopperSimulation(w0, w, th, g, b);
		this->v0PrintObject << this->kl << "," << this->gam << "," << g << "," << w_scale << "," \
			<< result << "," << this->kb << "," << this->ka << "," << this->spK << "," << start << \
			"," << interval << "," << this->coefv << "," << this->coefu << endl;
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
#if optimizer
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
	opt.set_xtol_rel(1e-4);
	std::vector<double> x(4);
	x[0] = 1; x[1] = 0.1; x[2] = 0.2; x[3] = 0.4;
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
#endif