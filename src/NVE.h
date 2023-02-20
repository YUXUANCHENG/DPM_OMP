#ifndef NVE_H
#define NVE_H

#include "cellPacking2D.h"
#include <random>
#include <stack>
#include "deformableParticles2D.h"

extern bool replaceFlag;
extern bool variableExtFflag;
extern bool settleDown;
extern bool frictionFlag;

extern bool softFlag;

class DPMNVEsimulator{
public:
	cellPacking2D * cellpointer = nullptr;
	double init_E = 0, init_U = 0, init_K = 0;
	int dof = 2;

	DPMNVEsimulator(cellPacking2D* cell) {
		cellpointer = cell;
		dof = 2;
	}

	DPMNVEsimulator() = default;

	virtual void NVERoutine() {

		for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
			cellpointer->cell(ci).verletPositionUpdate(cellpointer->dt);
			cellpointer->cell(ci).updateCPos();
		}

		// reset contacts before force calculation
		cellpointer->resetContacts();
		// calculate forces
		cellpointer->calculateForces();
		// update velocities
		verletVelocityUpdate();
	}

	virtual void verletVelocityUpdate(){
		for (int ci = 0; ci < cellpointer->NCELLS; ci++)
			cellpointer->cell(ci).verletVelocityUpdate(cellpointer->dt);
	}

	virtual void NVEsimulation(double T, double v0, double t_scale, int frames) {
		injectT(v0);
		NVEsimulationNoInjection(T, v0, t_scale, frames);
	}

	virtual void NVEsimulationNoInjection(double T, double v0, double t_scale, int frames) {
		int count = 0;
		int print_frequency = floor(T / (cellpointer->dt0 * t_scale * frames));

		for (double t = 0.0; t < T; t = t + cellpointer->dt0 * t_scale) {
			printRoutine(count, print_frequency, t);
			NVERoutine();
			count++;
		}
	}

	virtual void printRoutine(int count, int print_frequency, double t){
		cellpointer->printRoutine(count, print_frequency, t, init_E, init_U);
	}

	virtual void injectT(double v0) {
		int ci;
		int totalNV = 0;
		double totalDof;
		for (ci = 0; ci < cellpointer->NCELLS; ci++) {
			//system_mass += cell(ci).getNV() * PI * pow(0.5 * cell(ci).getdel() * cell(ci).getl0(), 2);
			totalNV += cellpointer->cell(ci).getNV();
		}
		totalDof = 2 * dof * totalNV / cellpointer->NCELLS;
		double scaled_v = cellpointer->scale_v(v0);
		init_K = cellpointer->cal_temp(scaled_v);
		init_U = cellpointer->totalPotentialEnergy();
		init_E = totalDof * init_K + init_U;
		// Reset velocity
		cellpointer->resetV();
		cellpointer->rescal_V(init_E);
	}

};


class BumpyNVEsimulator : public DPMNVEsimulator {
public:

	BumpyNVEsimulator(cellPacking2D* cell){
		cellpointer = cell;
		dof = 3;
	}
	BumpyNVEsimulator() = default;

	virtual void NVERoutine() {
		// update positions
		cellpointer->spPosVerlet();
		cellpointer->bumpyRotation();
		// reset contacts before force calculation
		cellpointer->resetContacts();
		// calculate forces
		cellpointer->bumpy_Forces();
		// update velocities
		verletVelocityUpdate();
	}

	virtual void verletVelocityUpdate(){
		cellpointer->sp_VelVerlet();
		cellpointer->bumpy_angularV();
	}

	virtual void injectT(double v0) {
		double totalDof;
		totalDof = dof;
		double scaled_v = cellpointer->scale_v(v0);
		init_K = cellpointer->cal_temp(scaled_v);
		init_U = cellpointer->totalPotentialEnergy();
		init_E = totalDof * init_K + init_U;
		// Reset velocity
		cellpointer->resetV();
		cellpointer->rescal_V(init_E);
	}
};

class BumpyLangevin : public BumpyNVEsimulator {
public:
	// Mersenne twister PRNG, initialized with seed from previous random device instance
	std::mt19937 gen;
	std::normal_distribution<double> dist;

	BumpyLangevin(cellPacking2D* cell){
		cellpointer = cell;
		dof = 3;
		std::random_device rd;
		gen = std::mt19937(rd());
		dist = std::normal_distribution<double>(0, 1);
		
	}
	BumpyLangevin() = default;
	// virtual void NVEsimulationNoInjection(double T, double v0, double t_scale, int frames) {
	// 	cellpointer->dt /= 10;
	// 	BumpyNVEsimulator::NVEsimulationNoInjection(T, v0, t_scale, frames);
	// 	cellpointer->dt *= 10;
	// }
	virtual void verletVelocityUpdate(){
		cellpointer->sp_VelVerlet_Langevin(1e-2, 2* init_K/ cellpointer->NCELLS, dist, gen);
		cellpointer->bumpy_angularV();
	}
	// virtual void printRoutine(int count, int print_frequency, double t){
	// 	if (count % print_frequency == 0)
	// 		cout << "t = " << t << endl;
	// }
};

class DPMLangevin : public DPMNVEsimulator {
public:
	// Mersenne twister PRNG, initialized with seed from previous random device instance
	std::mt19937 gen;
	std::normal_distribution<double> dist;

	DPMLangevin(cellPacking2D* cell){
		cellpointer = cell;
		dof = 2;
		std::random_device rd;
		gen = std::mt19937(rd());
		dist = std::normal_distribution<double>(0, 1);
	}
	DPMLangevin() = default;

	virtual void verletVelocityUpdate(){
		for (int ci = 0; ci < cellpointer->NCELLS; ci++)
			cellpointer->cell(ci).velVerlet_Langevin(cellpointer->dt, 1e-2, 2* init_K/ cellpointer->NCELLS, dist, gen);
		cellpointer->conserve_momentum();
	}
	// virtual void printRoutine(int count, int print_frequency, double t){
	// 	if (count % print_frequency == 0)
	// 		cout << "t = " << t << endl;
	// }
};

class DPMhopperSimulator {
public:
	cellPacking2D* cellpointer = nullptr;
	int closed = 1;
	int N_inside = 0;
	int clog_count = 0;
	int flowCount = 0;
	int placementNumber = 0;
	bool startPullFlag = false;
	double w0, w, th, g, b;
	int t;
	int lastFlowTime = 0;
	std::ofstream intervalRobj;

	DPMhopperSimulator(cellPacking2D* cell) {
		cellpointer = cell;
	}

	DPMhopperSimulator() = default;

	int hopperFlow(double w0, double w, double th, double g, double b) {
		this->w0 = w0; this->w = w; this->th = th; this->g = g; this->b = b;
		std::ofstream flowRobj;
		flowRobj.open("flowRate.txt");
		intervalRobj.open("flowInterval.txt");
		int result;
		int flowFlag = 0;
		for (t = 0; t < round(cellpointer->NT * 0.005 / cellpointer->dt); t++) {
			bool startFlag;
			if (replaceFlag)
				startFlag = (t > cellpointer->NT / 500 || Ke() < 1e-4 * N_inside);
			else 
				startFlag = (t > 1e4 && Ke() < 1e-4 * N_inside);
				// startFlag = (t > cellpointer->NT / 500 && Ke() < 1e-4 * N_inside);

			if (closed == 1 && startFlag) 
				closed = 0;

			addBack();
			cellpointer->printRoutine(t, round(cellpointer->NPRINT * 0.005 / cellpointer->dt), t, N_inside, closed);
			if (replaceFlag && closed == 0 && (t+1) % int (1e4 * 0.005 / cellpointer->dt) == 0) {
			// if (replaceFlag && closed == 0 && (t+1) % (cellpointer->NPRINT*1) == 0) {
				// addBack();
				if (flowCount)
					flowFlag = 0;
				else 
					flowFlag++;
				double rate = (double) flowCount / 1e4;
				// double rate = (double) flowCount / (double) cellpointer->NPRINT;
				flowRobj << rate << endl;
				flowCount = 0;
			}

			if (variableExtFflag)
			{
				if (t < 100000)
				{
					g = 0;
				}
				else if (t%1000 == 0)
				{
					startPullFlag = true;
					g += 0.01;
					// flowRobj << g << "," << cellpointer->cell(0).cforce(0)  << endl;
					flowRobj << g << "," << Ke()  << endl;
				}
			}

			hopperRoutine(w0, w, th, g, b);
			result = checkTermination();
			// if (result < 2 || flowFlag > 6)
			if (result < 2 )
			{	
				cellpointer->printRoutine(0, round(cellpointer->NPRINT * 0.005 / cellpointer->dt), t, N_inside, closed);
				return result;
			}
			// if (flowFlag > 6)
			// {	
			// 	cellpointer->printRoutine(0, round(cellpointer->NPRINT * 0.005 / cellpointer->dt), t, N_inside, closed);
			// 	return 0;
			// }
		}
		cout << "time out" << endl;
		cellpointer->printRoutine(0, round(cellpointer->NPRINT * 0.005 / cellpointer->dt), t, N_inside, closed);
		return 1;
	}

	int checkTermination() {
		if (N_inside == 0 || cellpointer->pistonX > cellpointer->L.at(0))
		{
			return 0;
		}
		if (cellpointer->NCELLS == 1 && Ke() > 0.1 && startPullFlag)
			return 0;

		// if (cellpointer->NCELLS > 1 && Ke() < 1e-12 * N_inside * pow(cellpointer->cell(0).NV / 16.0, 2) && closed == 0)
		if (Ke() < 1e-12 * N_inside * pow(cellpointer->cell(0).NV / 16.0, 2) && closed == 0)
		{
			clog_count++;
			if (clog_count > 100)
				return 1;
		}
		else
			clog_count = 0;
		return 2;
	}

	virtual double Ke() {
		return cellpointer->totalKineticEnergy();
	}

	virtual void hopperRoutine(double w0, double w, double th, double g, double b) {
		N_inside = 0;
		if (softFlag){
			cellpointer->hopperPosVerletSP();
		}
		for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
			if (cellpointer->cell(ci).inside_hopper)
			{
				N_inside++;
				// cellpointer->cell(0).inside_hopper = 0;
				if (!softFlag)
				{
					cellpointer->cell(ci).verletPositionUpdate(cellpointer->dt);
					cellpointer->cell(ci).updateCPos();
					if (ci < Nobs)
					{
						// reset
						deformableParticles2D * currentCell = &(cellpointer->cell(ci));
						for (int vi=0; vi<currentCell->getNV(); vi++){
							for (int d=0; d<cellpointer->NDIM; d++)
								currentCell->setVPos(vi,d,obstacleArray.at(ci).at(d) + currentCell->vrel(vi,d));
						}
						cellpointer->cell(ci).updateCPos();
					}
				}
				// if still inside hopper
				// int checkEk = cellpointer->cell(ci).totalKineticEnergy() > 1e3;
				// if (checkEk)
				// {
				// 	cout << "Ek error" << endl;
				// 	cout << ci << endl;
				// 	cellpointer->printRoutine(0, cellpointer->NPRINT, 0, N_inside, closed);
				// 	exit(0);
				// }
				if (cellpointer->cell(ci).cpos(0) > cellpointer->L.at(0) * 1.4 && cellpointer->cell(ci).cpos(0) > cellpointer->L.at(0) + 5 * sqrt(cellpointer->cell(ci).geta0()/PI))
				//if (cellpointer->cell(ci).cpos(0) > cellpointer->L.at(0) * 1.5 || cellpointer->cell(ci).cpos(0) < cellpointer->BoundaryCoor.at(0))
				{
					cellpointer->cell(ci).inside_hopper = 0;
					cellpointer->cell(ci).setCVel(0,0);
					cellpointer->cell(ci).setCVel(1,0);
					cellpointer->cell(ci).setCPos(0,100);
					cellpointer->cell(ci).setCPos(1,cellpointer->L.at(1)/2);
					cellpointer->cell(ci).regularPolygon();
				}
			}
		}
		// reset contacts before force calculation
		cellpointer->resetContacts();
		// calculate forces
		if (softFlag)
			cellpointer->hopperForcesSP(cellpointer->diskRadii,w0,w,th,g,closed);
		else
			cellpointer->hopperForces(w0, w, th, g, closed);
		// update velocities
		if (softFlag)
			cellpointer->hopperVelVerletSP(cellpointer->diskRadii, b);
		else
		{
			for (int ci = 0; ci < cellpointer->NCELLS; ci++)
				if (cellpointer->cell(ci).inside_hopper)
				{
					cellpointer->cell(ci).verletVelocityUpdate(cellpointer->dt, b);
				}
		}
	}
	double getMaxHight(double y)
	{
		double maxHight = 10;
		int maxIndex = 0;
		for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
			if (cellpointer->cell(ci).cpos(1) < y + 1.5 && cellpointer->cell(ci).cpos(1) > y - 1.5 ){
				if (cellpointer->cell(ci).cpos(0) < maxHight)
				{
					maxHight = cellpointer->cell(ci).cpos(0);
					maxIndex = ci;
				}
			}
		}
		for (int vi = 0; vi < cellpointer->cell(maxIndex).NV; vi ++)
		{
			double vY = cellpointer->cell(maxIndex).vpos(vi, 0);
			if ( vY < maxHight)
				maxHight = vY;
		}
		return maxHight;
	}
	double getVatMaxHeight(int cindex)
	{
		double maxHight = 0;
		int maxIndex = -1;
		double y = cellpointer->cell(cindex).cpos(1);
		double x = cellpointer->cell(cindex).cpos(0);
		for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
			if (cellpointer->cell(ci).inside_hopper == 1 && cellpointer->cell(ci).cpos(1) < y + 1.5 && cellpointer->cell(ci).cpos(1) > y - 1.5 ){
				if (cellpointer->cell(ci).cpos(0) < maxHight && cellpointer->cell(ci).cpos(0) > x)
				{
					maxHight = cellpointer->cell(ci).cpos(0);
					maxIndex = ci;
				}
			}
		}
		return maxIndex >= 0? cellpointer->cell(maxIndex).cvel(0) : 0;
	}
	void addBack()
	{
		int outside = 0;
		int inBetween = 1;
		double meanV = 0;
		double maxV = -10;
		stack<deformableParticles2D *> stack;
		for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
			if (cellpointer->cell(ci).inside_hopper == 0)
			{
				outside ++;
				flowCount ++;
				stack.push(&(cellpointer->cell(ci)));
			}
			else
			{
				if (cellpointer->cell(ci).cvel(0) > maxV)
					maxV = cellpointer->cell(ci).cvel(0);
				meanV += cellpointer->cell(ci).cvel(0);
				// if (cellpointer->cell(ci).cpos(0) > -21 && cellpointer->cell(ci).cpos(0) < -18)
				// {
				// 	meanV += cellpointer->cell(ci).cvel(0);
				// 	inBetween ++;
				// }
			}
		}
		
		if (outside > 0)
		{
			int interval = t - lastFlowTime;
			lastFlowTime = t;
			intervalRobj << interval << endl;

			meanV /= (cellpointer->NCELLS - outside);
			// meanV /= inBetween;
			int safeZone = 2;
			// int safeZone = 15;
			int NperLine = floor(cellpointer->L.at(1)/2) - safeZone;
			while (!stack.empty())
			{
				for(int i = 0; i < NperLine; i++)
				{
					placementNumber = placementNumber%NperLine;
					double displace = round((placementNumber+1e-4)/2.0) * 2 * pow(-1,placementNumber%2);

					if (!stack.empty())
					{
						placementNumber ++;
						deformableParticles2D * currentCell = stack.top();
						double maxHight = getMaxHight(displace + cellpointer->L.at(1)/2);
						currentCell->setCPos(0,maxHight - 1);
						currentCell->setCPos(1,displace + cellpointer->L.at(1)/2);
						currentCell->regularPolygon();
						// if (frictionFlag)
							currentCell->setCVel(0,meanV);
						// else
							// currentCell->setCVel(0,maxV);
						currentCell->setCVel(1,0);

						// update real-space positions
						for (int vi=0; vi<currentCell->getNV(); vi++){
							for (int d=0; d<cellpointer->NDIM; d++)
								currentCell->setVPos(vi,d,currentCell->cpos(d) + currentCell->vrel(vi,d));
						}

						currentCell->inside_hopper = 2;
						stack.pop();
					}
				}
			}
			settleDown = 1;
			//settleParticles();
			settleDown = 0;
			for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
				if (cellpointer->cell(ci).inside_hopper == 2)
				{
					double v = getVatMaxHeight(ci);
					double vFactor = 1;
					double bThreash = 5e-3;
					if (this->b < bThreash && (v*v-2*g*1.5/16)>0)
					{
						vFactor = sqrt(v*v-2*g*1.5/16)/v;
						// vFactor = meanV/v;
					}

					cellpointer->cell(ci).setCVel(0,v*vFactor);
					cellpointer->cell(ci).setCVel(1, 0.1*v*pow(-1,ci));
					cellpointer->cell(ci).inside_hopper = 1;
				}
			}
		}
	}
	void settleParticles(){
		bool result = false;
		double maxT, Tlim;
		// if (frictionFlag)
		// 	maxT = 1e4;
		// else
		// 	maxT = 1e3;
		if (b < 0.01)
			maxT = 1e4;
		else
			maxT = 1e3 * 0.1/b;
		Tlim = maxT;
		// if (t>2e5)
		// 	Tlim = 3 * maxT;
		for (int time = 0; time < Tlim; time++) {
			
			for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
				if (cellpointer->cell(ci).inside_hopper == 2)
				{
					cellpointer->cell(ci).verletPositionUpdate(cellpointer->dt);
					cellpointer->cell(ci).updateCPos();
				}
			}
			// reset contacts before force calculation
			cellpointer->resetContacts();
			// calculate forces
			cellpointer->hopperForces(w0, w, th, g, closed);
			// update velocities
			for (int ci = 0; ci < cellpointer->NCELLS; ci++)
				if (cellpointer->cell(ci).inside_hopper == 2)
					cellpointer->cell(ci).verletVelocityUpdate(cellpointer->dt, 0.01);
			
			int nSettling = 0;
			
			if (time%1000 == 0)
			{
				double tailEk = 0;
				for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
					if (cellpointer->cell(ci).inside_hopper == 2)
					{
						nSettling += 1;
						tailEk += cellpointer->cell(ci).totalKineticEnergy();
					}
				}
				result = tailEk < 1e-5 * nSettling;
				cout << "settling particles Ek = " << tailEk << endl;
			}
			if (result)
				return;
		}
	}
};
class BumpyHopperSimulator : public DPMhopperSimulator {
public:

	BumpyHopperSimulator(cellPacking2D* cell) {
		cellpointer = cell;
	}

	BumpyHopperSimulator() = default;

	virtual double Ke() {
		return cellpointer->totalKineticEnergy() + cellpointer->totalRotaionalK();
	}
	virtual void hopperRoutine(double w0, double w, double th, double g, double b) {
		N_inside = 0;
		for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
			if (cellpointer->cell(ci).inside_hopper)
				N_inside++;
		}
		cellpointer->spPosVerlet();
		cellpointer->bumpyRotation();
		for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
			if (cellpointer->cell(ci).cpos(0) > cellpointer->L.at(0) * 1.4 && cellpointer->cell(ci).cpos(0) > cellpointer->L.at(0) + 5 * sqrt(cellpointer->cell(ci).geta0()/PI))
				cellpointer->cell(ci).inside_hopper = 0;
		}
		//addBack();
		// reset contacts before force calculation
		cellpointer->resetContacts();
		// calculate forces
		cellpointer->hopperForces(w0, w, th, g, closed);
		// update velocities
		cellpointer->sp_VelVerlet(b);
		cellpointer->bumpy_angularV(b);

	}
};

#endif
