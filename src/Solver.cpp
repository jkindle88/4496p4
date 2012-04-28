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
	int MAX_FRAMES = 100;

	// loop over the frames
	for (int f = 0; f < MAX_FRAMES; f++)
	{
		cout << "starting frame: " << f << endl;

		// calculate constraints for the given frame
		computeConstraints(f,true);

		// compute objective function
		double objectiveFunction = EvaluateObjectiveFunction(frameCounter);

		double localStepSize = mStepSize;
		double localEpsilon = mEpsilon;
		int iterations = 0;

		int decreaseStepCounter = 0;

		while (objectiveFunction > localEpsilon )//&& iterations < mMaxIterations)
		{
			std::cout << "Iteration " << iterations << "\tObjective Function: " << objectiveFunction 
						<< "\tStep Size: " << localStepSize << std::endl;


			// calculate gradient
			Vecd gradient = computeGradient(f);

			// get old dofs
			Vecd oldDofs;
			oldDofs.SetSize(mModel->GetDofCount());
			mModel->mDofList.GetDofs(&oldDofs);

			// move dofs
			Vecd newDofs = oldDofs - localStepSize * gradient;

			// update dofs
			mModel->SetDofs(newDofs);

			// recompute constraint values
			computeConstraints(f,false);
			// calculate new objective function value
			double newObjectiveFunction = EvaluateObjectiveFunction(frameCounter);
			// make sure we always decrease so that we don't get stuck in some infinite loop
			if (newObjectiveFunction < objectiveFunction)
			{
				// adaptive step size
				// if difference between objective functions is greater than epsilon
				// and certain number of iterations have passed, increase step size
				if (newObjectiveFunction > mEpsilon &&
					iterations > 0 && iterations % 20 == 0)
				{
					localStepSize *= 4.0;
				}

				objectiveFunction = newObjectiveFunction;
			}
			else
			{

				// half step size
				decreaseStepCounter++;
				localStepSize /= 2.0;

				//localEpsilon *= 2.0;
				std::cout << "OBJECTIVE FUNCTION INCREASED!\t" << newObjectiveFunction << "\tEpsilon " << localEpsilon << std::endl;
				//if (decreaseStepCounter > 10)	// if we reach this pont, just go ahead and stop
				//{
					//objectiveFunction = 0.0;
				//}
				//else
				//{
					if (decreaseStepCounter < 300)
					{
						// reset to old dofs
						mModel->SetDofs(oldDofs);
					}
					else
					{
						// at this point, give up; it isn't worth it
						objectiveFunction = 0.0;
					}
				//}
			}

			// update iteration counter
			iterations++;
		}
		
		//UI->mGLWindow->flush();	// update screen
		//std::cout << "Ending frame " << frameCounter << std::endl;

	}
}

/*
* f - the current frame counter
* create - if true, create the constraints
*/
void Solver::computeConstraints(int f,bool create) {

	//Get the list of marker handles on the model
	Markers &handles = mModel->mHandleList;
	
	for (int i = 0; i < handles.size(); i++)
	{
		// get handle on model and contraint position
		Marker *handle = handles[i];
		Vec3d &c = mModel->mOpenedC3dFile->GetMarkerPos(f, i);

		// if constraint is (0, 0, 0), then we can disregard it
		if ( !(c[0] == 0 && c[1] == 0 && c[2] == 0) )
		{
			mConstraintList.push_back(Constraint(mModel, i));
		}
	}


}


