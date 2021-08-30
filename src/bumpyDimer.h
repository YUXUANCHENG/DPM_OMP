#ifndef BUMPYDIMER_H
#define BUMPYDIMER_H

#include "bumpyEllipse.h"
#include "NVE.h"
#include "cellPacking2D.h"

class BumpyDimer : public virtual BumpyEllipse {
public:

	//using BumpyEllipse::BumpyEllipse;
	BumpyDimer(int ncells, int nt, int nprint, double l, double s) :cellPacking2D::cellPacking2D(ncells, nt, nprint, l, s) {};
	BumpyDimer() = default;
	virtual void initializeGel(int NV, double phiDisk, double sizeDispersion, double delval) {
		// local variables
		int ci, vi, d, nvtmp;
		int lx, ly;
		double calA0tmp;
		double xpos, ypos, dx, dy;
		double xmin, xmax, ymin, ymax;
		double r1, r2, g1, areaSum;
		double rtmp, a0tmp, l0tmp, p0tmp;
		double angle = acos(ratio_a_b - 1);
		// seed random number generator
		////drand48()(23562457*seed);

		// minimum number of vertices
		const int nvmin = 12;

		// output to console
		cout << "		-- In gelation initialization, initializing cells and relaxing initial overlaps as repulsive SP particles" << endl;

		// initialize length scales as gaussian random variables (will becomes area square roots)
		areaSum = 0.0;
		vector<double> lenscales(NCELLS, 0.0);
		vector<double> diskradii(NCELLS, 0.0);
		for (ci = 0; ci < NCELLS; ci++) {
			// generate random numbers
			r1 = (double)rand() / (RAND_MAX + 1.0);
			r2 = (double)rand() / (RAND_MAX + 1.0);

			// calculate gaussian random variable using Box-Muller transform
			g1 = sqrt(-2.0 * log(r1)) * cos(2 * PI * r2);

			// get root area
			lenscales.at(ci) = g1 * sizeDispersion + 1.0;

			// effective particle area
			a0tmp = lenscales.at(ci) * lenscales.at(ci);

			// add to lenscales sum for boundary size
			areaSum += a0tmp;

			// save disk radius based on area square root
			diskradii.at(ci) = sqrt(ratio_a_b * a0tmp / PI);
		}

		// determine box length from particle sizes and input packing fraction
		for (d = 0; d < NDIM; d++)
			L.at(d) = sqrt(1.1 * areaSum * ratio_a_b / phiDisk);

		// set phi to input
		phi = phiDisk;

		// reseed rng
		////drand48()(56835698*seed);

		// initialize cell information
		cout << "		-- Ininitializing plant cell objects" << endl;
		for (ci = 0; ci < NCELLS; ci++) {
			// boundary information ( SET PBCS TO 1 for plant cells )
			for (d = 0; d < NDIM; d++) {
				cell(ci).setL(d, L.at(d));
				cell(ci).setpbc(d, 1);
			}

			// number of vertices ( SIGMA SETS # OF VERTS )
			nvtmp = NV;
			if (nvtmp > nvmin)
				cell(ci).setNV(nvtmp);
			else
				cell(ci).setNV(nvmin);

			// array information
			cell(ci).initializeVertices();
			cell(ci).initializeCell();

			// initialize cells as regular polygons
			calA0tmp = nvtmp * tan(PI / NV) / PI;

			// preferred area is lenscales squared
			a0tmp = lenscales.at(ci) * lenscales.at(ci);

			// initial length of polygon side
			double ration_square = ratio_a_b * ratio_a_b;
			l0tmp = ((4 * PI - 4 * angle) * diskradii.at(ci) / ratio_a_b) / nvtmp;

			// set preferred area and length 
			cell(ci).seta0(a0tmp);
			cell(ci).setl0(l0tmp);
			cell(ci).setdel(delval);
		}

		// initialize particle positions
		cout << "		-- Ininitializing cell positions" << endl;

		// set min and max values of positions
		xmin = 0;
		xmax = L.at(0);
		ymin = 0;
		ymax = L.at(1);

		// initialize positions of each cell
		for (ci = 0; ci < NCELLS; ci++) {
			// map onto random x and y position
			xpos = (xmax - xmin) * (double)rand() / (RAND_MAX + 1.0) + xmin;
			ypos = (ymax - ymin) * (double)rand() / (RAND_MAX + 1.0) + ymin;

			// set as initial position of com
			cell(ci).setCPos(0, xpos);
			cell(ci).setCPos(1, ypos);

			// initialize vertices as a regular polygon
			cell(ci).dimer(ratio_a_b);

			// perturb vertex positions a little bit
			//cell(ci).vertexPerturbation(0.1);
			cell(ci).printlengthscale(lengthscalePrintObject);
		}

		breakAlignment();

		// set time step
		vertexDPMTimeScale(0.1);

		// use FIRE in PBC box to relax overlaps
		cout << "		-- Using FIRE to relax overlaps..." << endl;
		for (ci = 0; ci < NCELLS; ci++)
			diskradii.at(ci) *= 1.1;
		fireMinimizeSP(diskradii);
		for (ci = 0; ci < NCELLS; ci++)
			diskradii.at(ci) /= 1.1;
		lengthscalePrintObject << L.at(0) << endl << L.at(1) << endl;

		// printJammedConfig_yc();
		// printCalA();
		// printContact();

	}

	virtual void initShape(double ratio, deformableParticles2D & cell){
		cell.dimer(ratio);
	}

	
};

// initialize vertex positions so cell begins as ellipses
void deformableParticles2D::dimer(double ratio) {
	// local variables
	int i, j;
	double angleArg = 0.0;
	double a = sqrt(ratio * a0 / PI);;
	double b = a / ratio;
	double angle = acos(ratio - 1);
	double k, x, y, sign;

	// check if NV has been set > 0
	if (NV <= 0) {
		cout << "	ERROR: in regularPolygon(), NV = " << NV << ", so not set properly. Ending." << endl;
		exit(1);
	}
	else if (a0 < 0.1) {
		cout << "	ERROR: in regularPolygon(), a0 = " << a0 << ", so too small and not set properly. Ending." << endl;
		exit(1);
	}

	// loop over vertices, set positions using rotations
	for (i = 0; i < int (NV / 2); i++) {
		angleArg = angle + i * 4 * ( PI - angle) / NV;
		setVRel(i, 0, -(a - b) + b * cos(angleArg));
		setVRel(i, 1, b * sin(angleArg));
	}
	for (i = int(NV / 2); i < NV; i++) {
		j = i - int(NV / 2);
		angleArg = PI + angle + j * 4 * (PI - angle) / NV;
		setVRel(i, 0, (a - b) + b * cos(angleArg));
		setVRel(i, 1, b * sin(angleArg));
	}
	// output
	cout << " 	-- creating regular polygon with a0 = " << a0 << ", area = " << polygonArea() << " and perimeter = " << perimeter() << ", so init calA0 = " << pow(perimeter(), 2.0) / (4.0 * PI * polygonArea()) << ", and calA = " << calA() <<", compare to " << NV * tan(PI / NV) / PI << endl;
}

#endif