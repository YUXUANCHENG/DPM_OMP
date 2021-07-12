#ifndef BUMPY_PARALLEL_H
#define BUMPY_PARALLEL_H

#include <chrono>
#include "DPM_Parallel.h"
#include "bumpy.h"
#include "NVE.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

class Bumpy_Parallel : public virtual Bumpy, public virtual DPM_Parallel {
public:
	// g++ does not support multiple inheritance of constructors???
	//using Bumpy::Bumpy;
	Bumpy_Parallel(int ncells, int nt, int nprint, double l, double s) :cellPacking2D::cellPacking2D(ncells, nt, nprint, l, s) {};
	Bumpy_Parallel() = default;

	using DPM_Parallel::initialize_subsystems;
	
	virtual void hopperForces(double w0, double w, double th, double g, int closed){
// auto t1 = high_resolution_clock::now();
		bumpy_Forces();
// auto t2 = high_resolution_clock::now();
		if (gOn){
#pragma omp parallel for
			for (int ci = 0; ci < NCELLS; ci++) 
				cell(ci).gravityForces(g, gDire);
		}
// auto t3 = high_resolution_clock::now();
		hopperWallForcesDP(w0,w,th,closed);
// auto t4 = high_resolution_clock::now();
// duration<double, std::milli> ms_double_1 = t2 - t1;
// duration<double, std::milli> ms_double_2 = t3 - t2;
// duration<double, std::milli> ms_double_3 = t4 - t3;
// std::cout << ms_double_1.count() << " " << ms_double_2.count() << " " << ms_double_3.count() << endl;
	}

	virtual void bumpy_Forces() { 
		// reset forces
		for (int ci = 0; ci < NCELLS; ci++) {
			// reset center of mass forces
			for (int d = 0; d < NDIM; d++)
			{
				cell(ci).setCForce(d, 0.0);
				cell(ci).torque = 0.0;
			}
			// reset vertex forces and interaction energy
			for (int vi = 0; vi < cell(ci).getNV(); vi++) {
				// forces
				for (int d = 0; d < NDIM; d++)
					cell(ci).setVForce(vi, d, 0.0);

				// energies
				cell(ci).setUInt(vi, 0.0);
			}
		}
		resetContacts();
// auto t1 = high_resolution_clock::now();
		subspaceManager();
// auto t2 = high_resolution_clock::now();
#pragma omp parallel for schedule (dynamic, 4)
		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
			subsystem[i].calculateForces_insub();
//auto t2 = high_resolution_clock::now();
#pragma omp parallel for schedule (dynamic, 4)
		for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
			subsystem[i].calculateForces_betweensub();
// auto t3 = high_resolution_clock::now();
// duration<double, std::milli> ms_double_1 = t2 - t1;
// duration<double, std::milli> ms_double_2 = t3 - t2;
// std::cout << ms_double_1.count() << " " << ms_double_2.count() << endl;

		for (int ci = 0; ci < NCELLS; ci++) {
			for (int d = 0; d < NDIM; d++)
				cell(ci).setCForce(d, cell(ci).cforce(d));
		}
		addUpStress();
	}

	virtual void printRoutine(int count, int print_frequency, double t, double init_E, double init_U)
	{
		Bumpy::printRoutine(count, print_frequency, t, init_E, init_U);
		if (count % print_frequency == 0)
			cout << "total vetex count = " << countTotalVertices() << endl;
	}

};

#endif