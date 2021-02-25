#include "bumpy.h"

void Bumpy::qsIsoCompression(double phiDisk, double deltaPhi, double Ftolerance) {
	int ci;
	double dr;
	for (ci = 0; ci < NCELLS; ci++) {
		cell(ci).setkint(1.0);
		cell(ci).setdel(1.0);
		cell(ci).seta(0.0);
		cell(ci).setCForce(0, 0.0);
		cell(ci).setCForce(1, 0.0);
		for (int vi = 0; vi < cell(ci).getNV(); vi++)
			cell(ci).setUInt(vi, 0.0);
	}

	//double phi = packingFraction();
	phi = packingFraction();

	if (phi > phiDisk){
		dr = sqrt((phi - 0.001) / phi);
		while (phi > phiDisk) {
			// scale lengths
			scaleLengths(dr);
			phi = packingFraction();
		}
	}
	else{
		// compute length scaler based on deltaPhi
		dr = sqrt((phi + 0.001) / phi);
		// loop until phi is the correct value
		while (phi < phiDisk) {
			// scale lengths
			scaleLengths(dr);
			phi = packingFraction();
			for (int ci = 0; ci < NCELLS; ci++)
				cell(ci).cal_inertia();
			// relax shapes (energies calculated in relax function)
			fireMinimize_bummpy();
		}
	}

	for (ci = 0; ci < NCELLS; ci++)
		cell(ci).cal_inertia();
	fireMinimize_bummpy();
}

void Bumpy::compressToInitial(double phiTarget, double deltaPhi, double Ftol) {
	qsIsoCompression(phiTarget, deltaPhi, Ftol);
	phi = packingFraction();
	if (phi > 0.6)
	{
		// loop until phi is the correct value
		for (int i = 0; i < 50; i++) {
			// scale lengths
			scaleLengths(0.99);
		}
		// loop until phi is the correct value
		for (int i = 0; i < 50; i++) {
			// scale lengths
			scaleLengths(1 / 0.99);
			for (int ci = 0; ci < NCELLS; ci++)
				cell(ci).cal_inertia();
			// relax shapes (energies calculated in relax function)
			fireMinimize_bummpy();
			phi = packingFraction();
		}
	}
}

void Bumpy::printRoutine(int count, int print_frequency, double t, double init_E, double init_U) {
	if (count % print_frequency == 0) {
		// calculate energies
		double U = totalPotentialEnergy();
		double K = totalKineticEnergy() + totalRotaionalK();
		//rescal_V(current_E);
		printJammedConfig_yc();
		phi = packingFraction();
		phiPrintObject << phi << endl;
		printCalA();
		printContact();
		printV();
		cout << "phi = " << phi << endl;
		cout << "E_INIT = " << init_E << " U_INIT = " << init_U << endl;
		cout << "E = " << U + K << " K = " << K << " Kr = " << totalRotaionalK() << endl;
		cout << "t = " << t << endl;
	}
}

void Bumpy::printV() {
	double v_x = 0;
	double v_y = 0;
	double energy = 0.0;
	double tot_k_energy = 0.0;

	for (int ci = 0; ci < NCELLS; ci++) {
		v_x = cell(ci).cal_mean_v(0);
		v_y = cell(ci).cal_mean_v(1);
		energy = 0.5 * (v_x*v_x + v_y*v_y) * cell(ci).getNV();
		tot_k_energy = cell(ci).totalKineticEnergy() + 0.5 * cell(ci).inertia * pow(cell(ci).angularV, 2);
		vPrintObject << v_x << "," << v_y <<  "," << energy/tot_k_energy << "," << tot_k_energy << endl;
	}
}

void Bumpy::resetV() {
	int ci, d;
	double rv;
	for (ci = 0; ci < NCELLS; ci++) {
		for (d = 0; d < NDIM; d++) {
			// get random direction
			rv = (double)rand() / (RAND_MAX + 1.0);
			cell(ci).setCVel(d, rv);
		}
	}
}

void Bumpy::NVEsimulation(double T, double v0, double t_scale, int frames) {
	LangevinSimulation(2000, v0, t_scale, 20);
	BumpyNVEsimulator simulator = BumpyNVEsimulator(this);
	simulator.NVEsimulationNoInjection(T, v0, t_scale, frames);
}

void Bumpy::LangevinSimulation(double T, double v0, double t_scale, int frames) {
	BumpyLangevin simulator = BumpyLangevin(this);
	simulator.NVEsimulation(T, v0, t_scale, frames);
}

void Bumpy::fireMinimizeF(double Ftol, double& Ftest, double& Ktest) {
	for (int ci = 0; ci < NCELLS; ci++)
		cell(ci).cal_inertia();
	fireMinimize_bummpy();
}