/**
* Solver.cpp
* @Author: Joe Kindle
* Worker class for solving inverse kinematics.
* Will also be used to store the relevant values.
*/

#include <float.h>
#include "Solver.h"
#include "Model.h"
#include "RealTimeIKui.h"

extern RealTimeIKUI *UI;

// Constructor
Solver::Solver(double eps, double step, int maxIterations, int maxFrames, Model *model)
{
	Solver::mEps = eps;
	Solver::mStep = step;
	Solver::mMaxIters = maxIterations;
	Solver::mMaxFrames = maxFrames;
	Solver::mModel = model;
}

// Destructor
Solver::~Solver()
{
	Solver::mEps = NULL;
	Solver::mStep = NULL;
	Solver::mMaxIters = NULL;
	Solver::mMaxFrames = NULL;
	Solver::mModel = NULL;
}

void Solver::solve() {
	
}


