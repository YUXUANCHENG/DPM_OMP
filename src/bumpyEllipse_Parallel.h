#ifndef BUMPYELLIPSE_PARALLEL_H
#define BUMPYELLIPSE_PARALLEL_H

#include "bumpy_Parallel.h"
#include "bumpyEllipse.h"

class BumpyEllipse_Parallel : public Bumpy_Parallel,  public BumpyEllipse{
public:
    //using Bumpy::Bumpy;
    BumpyEllipse_Parallel(int ncells, int nt, int nprint, double l, double s) :cellPacking2D::cellPacking2D(ncells, nt, nprint, l, s) {};
	BumpyEllipse_Parallel() = default;

};

#endif