#include "cellPacking2D.h"
#include "NVE.h"
#include "taoSolver.h"

void cellPacking2D::resetV() {
	int ci, vi, d;
	double rv;
	for (ci = 0; ci < NCELLS; ci++) {
		for (vi = 0; vi < cell(ci).getNV(); vi++)
		{
			for (d = 0; d < NDIM; d++) {
				// get random direction
				rv = (double)rand() / (RAND_MAX + 1.0);
				cell(ci).setVVel(vi, d, rv);
			}
		}
	}
}

void cellPacking2D::printRoutine(int count, int print_frequency, double t, double init_E, double init_U) {
	if (count % print_frequency == 0) {
		double U = totalPotentialEnergy();
		double K = totalKineticEnergy();		
		printSubRoutine(count, print_frequency);
		cout << "phi = " << phi << endl;
		cout << "E_INIT = " << init_E << " U_INIT = " << init_U << endl;
		cout << "E = " << U + K << " K = " << K << " U = " << U << endl;
		cout << "t = " << t << endl;
	}
}

void cellPacking2D::printSubRoutine(int count, int print_frequency)
{
	if (count % print_frequency == 0) {
		printJammedConfig_yc();
		phi = packingFraction();
		phiPrintObject << phi << endl;
		printCalA();
		printContact();
		printV();
	}
}

void cellPacking2D::NVEsimulation(double T, double v0, double t_scale, int frames) {
	LangevinSimulation(4000, v0, t_scale, 20);
	DPMNVEsimulator simulator = DPMNVEsimulator(this);
	simulator.NVEsimulationNoInjection(T, v0, t_scale, frames);
}

void cellPacking2D::LangevinSimulation(double T, double v0, double t_scale, int frames) {
	DPMLangevin simulator = DPMLangevin(this);
	simulator.NVEsimulation(T, v0, t_scale, frames);
}

void cellPacking2D::ActiveBrownianSimulation(double T, double v0, double Dr, double vtau, double t_scale, int frames)
{
	DPMActiveBrownian simulator = DPMActiveBrownian(this, Dr, vtau);
	simulator.NVEsimulation(T, v0, t_scale, frames);
}

int cellPacking2D::hopperSimulation(double w0, double w, double th, double g, double b) {
	DPMhopperSimulator simulator = DPMhopperSimulator(this);
	return simulator.hopperFlow(w0, w, th, g, b);
}

void cellPacking2D::compressToInitial(double phiTarget, double deltaPhi, double Ftol){
	//double phi = packingFraction();
	phi = packingFraction();
	double temp_phiTarget;
	if (phiTarget < 0.7)
		temp_phiTarget = phiTarget;
	else
		temp_phiTarget = 0.7 + 1e-10;
	// compute length scaler based on deltaPhi
	double dr = sqrt((phi - 0.01) / phi);
	// loop until phi is the correct value
	while (phi > temp_phiTarget) {
		// scale lengths
		scaleLengths(dr);
		phi = packingFraction();
	}

	qsIsoCompression(phiTarget, deltaPhi, Ftol);
}

double* cellPacking2D::NVE_tao(double T, double v0, double Dr, double vtau, double t_scale, int frames) {
	DPMNVEsimulator simulator = DPMNVEsimulator(this);
	TaoSolver taoSolver = TaoSolver(this, &simulator);
	return taoSolver.NVE_tao(T, v0, Dr, vtau, t_scale, frames);
}


double* cellPacking2D::NVT_tao(double T, double v0, double Dr, double vtau, double t_scale, int frames) {
	DPMActiveBrownian simulator = DPMActiveBrownian(this, Dr, vtau);
	simulator.injectT(v0);
	TaoSolver taoSolver = TaoSolver(this, &simulator);
	return taoSolver.NVE_tao(T, v0, Dr, vtau, t_scale, frames);
}



