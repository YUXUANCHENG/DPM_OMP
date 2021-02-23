#include "bumpy.h"

void Bumpy::qsIsoCompression(double phiDisk, double deltaPhi, double Ftolerance) {
	int ci;
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
	// compute length scaler based on deltaPhi
	double dr = sqrt((phi + 0.001) / phi);
	// loop until phi is the correct value
	while (phi < phiDisk) {
		// scale lengths
		scaleLengths(dr);
		phi = packingFraction();
	}
	dr = sqrt((phi - 0.001) / phi);
	while (phi > phiDisk) {
		// scale lengths
		scaleLengths(dr);
		phi = packingFraction();
	}

	if (phi > 0.6)
	{
		// loop until phi is the correct value
		for (int i = 0; i < 100; i++) {
			// scale lengths
			scaleLengths(0.999);
		}
		// loop until phi is the correct value
		for (int i = 0; i < 100; i++) {
			// scale lengths
			scaleLengths(1 / 0.999);
			for (ci = 0; ci < NCELLS; ci++)
				cell(ci).cal_inertia();
			// relax shapes (energies calculated in relax function)
			fireMinimize_bummpy();
			phi = packingFraction();
		}
	}

	for (ci = 0; ci < NCELLS; ci++)
		cell(ci).cal_inertia();
	fireMinimize_bummpy();
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
		//printCalA();
		printContact();
		printV();
		cout << "phi = " << phi << endl;
		cout << "E_INIT = " << init_E << " U_INIT = " << init_U << endl;
		cout << "E = " << U + K << " K = " << K << " Kr = " << totalRotaionalK() << endl;
		cout << "t = " << t << endl;
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