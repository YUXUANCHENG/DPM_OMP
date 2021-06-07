#ifndef BUMPYDIMER_PARALLEL_H
#define BUMPYDIMER_PARALLEL_H

#include "bumpy_Parallel.h"
#include "bumpyDimer.h"

class BumpyDimer_Parallel : public Bumpy_Parallel,  public BumpyDimer{
public:
    //using Bumpy::Bumpy;
    BumpyDimer_Parallel(int ncells, int nt, int nprint, double l, double s) :cellPacking2D::cellPacking2D(ncells, nt, nprint, l, s) {};
	BumpyDimer_Parallel() = default;

};

#endif