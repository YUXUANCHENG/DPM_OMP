#ifndef CLI_HOPPER_H
#define CLI_HOPPER_H

#include "CLI.h"
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
	double g = 0.05;
	vector<double> radii;
	double w_scale = 2;
	double w;
	double gamafactor1 = 0;
	double gamafactor2 = 0;
	bool deformFlag = false;

	DPM_Hopper_CLI() {
		this->NT = 5e6;			// number of time steps for flow simulation
		this->NPRINT = 1e3;			// number of steps between printing
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
		this->timeStepMag = 0.00005;	
		//this->timeStepMag = 0.0005;	
		this->kl = 1;
		this->ka = 1;
		this->kb = 0;
		this->gamafactor1 = 0.2;
		// this->gamafactor2 = -3.3;
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
		std::ofstream deformationPrint;
		deformationPrint.open("deformation.txt");
		deformationPrint << (originalHeight - endHeight)/originalHeight << "," << angle << "," << length <<endl;
		deformationPrint.close();

	}

	void findParameter(char const* argv[]){
		double presetCalA = 1.14;
		double presetAngle = 17.9;
		double threashold = 0.3;

		deformFlag = true;
		this->NT = 5e6;
		this->th = PI/4;
		//this->timeStepMag = 0.00002;	
		//this->timeStepMag = 0.00005;	
		this->timeStepMag = 0.0001;	
		for (this->kl = 1; this->kl < 3; this->kl += 0.5){
			for (this->gamafactor1 = 0.01; this->gamafactor1 < 0.2; this->gamafactor1 += 0.05){
				for (this->gamafactor2 = -0.01; this->gamafactor2 > -0.2; this->gamafactor2 -= 0.05){
					this->ka = 1;
					this->kb = 0;
					
					// this->gamafactor2 = -3.3;
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
					if (this->particles->matchPreset(presetCalA, presetAngle, threashold)){
						cout << "Parameters found!" << endl;
						cout << "kl = " << this->kl << ", gama1 = " << this->gamafactor1 << ", gama2 = " << this->gamafactor2 << endl;
						return;
					}
				}
			}
		}
		cout<< "Parameters for preset shape were not found" << endl;
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

#endif