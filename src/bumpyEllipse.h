#ifndef BUMPYELLIPSES_H
#define BUMPYELLIPSES_H

#include "bumpy.h"
#include "NVE.h"
#include "DPM_Parallel.h"
#include "deformableParticles2D.h"

class BumpyEllipse : public virtual Bumpy {
public:
	double ratio_a_b = 1.8;

	//using Bumpy::Bumpy;
    BumpyEllipse(int ncells, int nt, int nprint, double l, double s) :cellPacking2D::cellPacking2D(ncells, nt, nprint, l, s) {};
	BumpyEllipse() = default;
	virtual void setRatio(double calA0) {
		double xpos = (double)rand() / (RAND_MAX + 1.0);
		deformableParticles2D celltemp = deformableParticles2D();
		// set as initial position of com
		celltemp.setNV(16);
		for (int d = 0; d < NDIM; d++) {
			celltemp.setL(d, 5);
			celltemp.setpbc(d, 1);
		}
		// array information
		celltemp.initializeVertices();
		celltemp.initializeCell();
		celltemp.seta0(1);
		celltemp.setCPos(0, xpos);
		celltemp.setCPos(1, xpos);
		// initialize vertices as a regular polygon
		initShape(ratio_a_b, celltemp);
		double currentCalA = celltemp.calA();
		while (!(currentCalA > calA0 * 0.999 && currentCalA < calA0 * 1.001))
		{
			if (currentCalA > calA0)		
				ratio_a_b *= 0.999;
			else
				ratio_a_b *= 1.001;
			initShape(ratio_a_b, celltemp);
			currentCalA = celltemp.calA();
		}
	}
	virtual void initializeGel(int NV, double phiDisk, double sizeDispersion, double delval) {
		// local variables
		int ci, vi, d, nvtmp;
		int lx, ly;
		double calA0tmp;
		double xpos, ypos, dx, dy;
		double xmin, xmax, ymin, ymax;
		double r1, r2, g1, areaSum;
		double rtmp, a0tmp, l0tmp, p0tmp;

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
			if (ci < NCELLS / 2)
				lenscales.at(ci) *= 1.4;

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
			nvtmp = ceil(lenscales.at(ci) * NV);
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
			double ration_square = ratio_a_b* ratio_a_b;
			l0tmp = (2 * PI * diskradii.at(ci) * sqrt((1 + ration_square)/ (2 * ration_square))) / nvtmp;

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
			cell(ci).ellipse(ratio_a_b);

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
		cell.ellipse(ratio);
	}

	void breakAlignment() {
		double theta;
		double xtemp, ytemp;

		for (int ci = 0; ci < NCELLS; ci++) {
			// update new position based on acceleration
			theta = 2* PI * (double)rand() / (RAND_MAX + 1.0);

			// rotate vertex
			for (int vi = 0; vi < cell(ci).getNV(); vi++)
			{
				xtemp = cell(ci).cpos(0) + cell(ci).vrel(vi, 0) * cos(theta) - cell(ci).vrel(vi, 1) * sin(theta);
				ytemp = cell(ci).cpos(1) + cell(ci).vrel(vi, 0) * sin(theta) + cell(ci).vrel(vi, 1) * cos(theta);
				cell(ci).setVPos(vi, 0, xtemp);
				cell(ci).setVPos(vi, 1, ytemp);
			}
		}
	}
};
// https://stackoverflow.com/questions/6972331/how-can-i-generate-a-set-of-points-evenly-distributed-along-the-perimeter-of-an
double ComputeArcOverAngle(double r1, double r2, double angle, double angleSeg);

double GetLengthOfEllipse(double deltaAngle, double a, double b)
{
	// Distance in radians between angles
	double numIntegrals = round(PI * 2.0 / deltaAngle);
	double length = 0;

	// integrate over the elipse to get the circumference
	for (int i = 0; i < numIntegrals; i++)
	{
		length += ComputeArcOverAngle(a, b, i * deltaAngle, deltaAngle);
	}

	return length;
}

double GetAngleForArcLengthRecursively(double currentArcPos, double goalArcPos, double angle, double angleSeg, double a, double b)
{
	double ARC_ACCURACY = 0.1;
	// Calculate arc length at new angle
	double nextSegLength = ComputeArcOverAngle(a, b, angle + angleSeg, angleSeg);

	// If we've overshot, reduce the delta angle and try again
	if (currentArcPos + nextSegLength > goalArcPos) {
		return GetAngleForArcLengthRecursively(currentArcPos, goalArcPos, angle, angleSeg / 2, a, b);

		// We're below the our goal value but not in range (
	} else if (currentArcPos + nextSegLength < goalArcPos - ((goalArcPos - currentArcPos) * ARC_ACCURACY)) {
		return GetAngleForArcLengthRecursively(currentArcPos + nextSegLength, goalArcPos, angle + angleSeg, angleSeg, a, b);

		// current arc length is in range (within error), so return the angle
	} else
		return angle;
}

double ComputeArcOverAngle(double r1, double r2, double angle, double angleSeg)
{
	double distance = 0.0;

	double dpt_sin = pow(r1 * sin(angle), 2.0);
	double dpt_cos = pow(r2 * cos(angle), 2.0);
	distance = sqrt(dpt_sin + dpt_cos);

	// Scale the value of distance
	return distance * angleSeg;
}

// initialize vertex positions so cell begins as ellipses
// void deformableParticles2D::ellipse(double ratio) {
// 	// local variables
// 	int i;
// 	double angleArg = 0.0;
// 	double a = sqrt(ratio * a0 / PI);;
// 	double b = a / ratio;
// 	double k, x, y, sign;

// 	// check if NV has been set > 0
// 	if (NV <= 0) {
// 		cout << "	ERROR: in regularPolygon(), NV = " << NV << ", so not set properly. Ending." << endl;
// 		exit(1);
// 	}
// 	else if (a0 < 0.1) {
// 		cout << "	ERROR: in regularPolygon(), a0 = " << a0 << ", so too small and not set properly. Ending." << endl;
// 		exit(1);
// 	}

// 	// loop over vertices, set positions using rotations
// 	for (i = 0; i < NV; i++) {
// 		angleArg = (2.0 * PI * i) / NV;
// 		k = tan(angleArg);
// 		sign = cos(angleArg) / (abs(cos(angleArg)) + 1e-30);
// 		x = sign * a / sqrt(1 + ratio * ratio * k * k);
// 		y = k * x;
// 		setVRel(i, 0, x);
// 		setVRel(i, 1, y);
// 	}
// 	// output
// 	cout << " 	-- creating regular polygon with a0 = " << a0 << ", area = " << polygonArea() << " and perimeter = " << perimeter() << ", so init calA0 = " << pow(perimeter(), 2.0) / (4.0 * PI * polygonArea()) << ", and calA = " << calA() <<", compare to " << NV * tan(PI / NV) / PI << endl;
// }

void deformableParticles2D::ellipse(double ratio)
{
	double a = sqrt(ratio * a0 / PI);;
	double b = a / ratio;

	// Distance in radians between angles measured on the ellipse
	double deltaAngle = 0.001;
	double circumference = GetLengthOfEllipse(deltaAngle, a, b);

	double arcLength = circumference/ NV;

	double angle = 0;

	// Loop until we get all the points out of the ellipse
	for (int numPoints = 0; numPoints < NV; numPoints++)
	{
		angle = GetAngleForArcLengthRecursively(0, arcLength, angle, deltaAngle, a, b);

		double x = a * cos(angle);
		double y = b * sin(angle);
		setVRel(numPoints, 0, x);
		setVRel(numPoints, 1, y);
	}
	cout << " 	-- creating regular polygon with a0 = " << a0 << ", area = " << polygonArea() << " and perimeter = " << perimeter() << ", so init calA0 = " << pow(perimeter(), 2.0) / (4.0 * PI * polygonArea()) << ", and calA = " << calA() <<", compare to " << NV * tan(PI / NV) / PI << endl;

}



#endif