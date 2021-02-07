#include "bumpy.h"

void Bumpy::initializeGel(int NV, double phiDisk, double sizeDispersion, double delval) {

	cellPacking2D::initializeGel(NV, phiDisk, sizeDispersion, delval);
	int ci;
	double phi = packingFraction();
	// compute length scaler based on deltaPhi
	double dr = sqrt((phi - 0.01) / phi);
	// loop until phi is the correct value
	while (phi > phiDisk) {
		// scale lengths
		scaleLengths(dr);
		phi = packingFraction();
	}
	if (phi > 0.85)
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
		}
	}

	for (ci = 0; ci < NCELLS; ci++)
		cell(ci).cal_inertia();
	fireMinimize_bummpy();

	for (ci = 0; ci < NCELLS; ci++) {
		cell(ci).setkint(1.0);
		cell(ci).setdel(1.0);
		cell(ci).seta(0.0);
		cell(ci).setCForce(0, 0.0);
		cell(ci).setCForce(1, 0.0);
		for (int vi = 0; vi < cell(ci).getNV(); vi++)
			cell(ci).setUInt(vi, 0.0);
	}
	bumpy_Forces();
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
	BumpyNVEsimulator simulator = BumpyNVEsimulator(this);
	simulator.NVEsimulation(T, v0, t_scale, frames);
}

void Bumpy::fireMinimizeF() {
	fireMinimize_bummpy();
}