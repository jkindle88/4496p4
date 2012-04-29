/**
* Solver.cpp
* @Author: Joe Kindle
* Worker class for solving inverse kinematics.
* Will also be used to store the relevant values.
*/

#include <float.h>
#include <map>
#include "Solver.h"
#include "Model.h"
#include "Marker.h"
#include "RealTimeIKui.h"		
#include "PhylterGLWindow.h"
#include "C3dFileInfo.h"
#include "Transform.h"
#include "TransformNode.h"

extern RealTimeIKUI *UI;	//Gives us access to UI->mGLWindow->flush()

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
	// Utility constants
	const int sIFreq = 20;		//STEP_INCREASE_FREQUENCY
	const double sIF = 4.0;	//STEP_INCREASE_FACTOR
	const double sDF = 2.0;	//STEP_DECREASE_FACTOR
	const int eIFreq = 2;		//EPSILON_INCREASE_FREQUENCY
	const double eIF = 2.0;	//EPSILON_INCREASE_FACTOR

	// Loop over all valid frames
	for (int f = 0; f < mMaxFrames; f++)
	{
		// Build C[]
		buildConstraintMat(f);

		//Objective function value
		double o = 0.0;

		Markers &handles = mModel->mHandleList;
		for (int i = 0; i < constraintIDs.size(); i++)
		{
			// Get handle and target position for the given constraint
			Marker *h = handles[constraintIDs[i]];
			Vec3d cPos = mModel->mOpenedC3dFile->GetMarkerPos(f, constraintIDs[i]);

			// Calculate values and store
			Vec3d cVal = h->mGlobalPos - cPos;
			double cLength = sqrlen(cVal);

			//Store results
			constraintVals.push_back(cVal);
			constraintLengths.push_back(cLength);

			//Update objective function
			o += cLength; 
		}
		
		// Main loop
		double e = mEps;
		double s = mStep;
		int iterCount = 0;
		int stepCount = 0;
		while (o > e) //while Objective Function > Epsilon
		{
			// Compute gradient -- TODO: Move this somewhere else for the sake of keeping this compact!
			Vecd gradient;
			gradient.SetSize(mModel->GetDofCount());
			gradient.MakeZero();

			for (int i = 0; i < constraintIDs.size(); i++)
			{
				//Get Constraint
				Vec4d cVec(constraintVals[i], 1.0);

				//Get Jacobian
				Matd J = computeJ(mModel->mHandleList, constraintIDs[i]);

			}

			// Get previous DOFs
			Vecd prevDofs;
			prevDofs.SetSize(mModel->GetDofCount());
			mModel->mDofList.GetDofs(&prevDofs);

			// Update DOFs
			Vecd nextDofs = prevDofs - s * gradient;
			mModel->SetDofs(nextDofs);

			// Recalculate constraints
			CalculateConstraints(f);
			// Recalculate objective function
			double newObjectiveFunction = EvaluateObjectiveFunction(f);

			/*
			if (newObjectiveFunction < objectiveFunction)
			{
				
				if (newObjectiveFunction > mEps && iterCount > 0 && iterCount % sIFreq == 0)
				{
					s *= sIF;
				}

				objectiveFunction = newObjectiveFunction;
				
			}
			else
			{
				stepCount++;

				// Decrease step size
				s = s / sDF;
				
				// based on how many times we've had to decrease the step, choose different options
				if (stepCount < mMaxIterations)
				{
					// Reset to the previous dofs and try again
					mModel->SetDofs(prevDofs);
				}
				else if (stepCount < eIFreq*mMaxIterations)
				{
					// Try increasing epsilon
					localEpsilon *= mEpsilonIncreaseFactor;

					// Reset to the previous dofs and try again
					mModel->SetDofs(prevDofs);
				}
				else
				{
					// Otherwise, give up
					objectiveFunction = 0.0;
				}
			}*/

			// Update iteration counter
			iterCount++;
		}

		// Update the frame counter
		UI->mFrameCounter_cou->value(f);

		//Refresh the screen
		UI->mGLWindow->flush();
	}	
}

/**
* Builds the Constraint Matrix C[]
*/
void Solver::buildConstraintMat(int f) {
	// Get the handles on the model and loop
	Markers &handles = mModel->mHandleList;
	for (int i = 0; i < handles.size(); i++)
	{
		// Evaluate the constraint
		Marker *h = handles[i];
		Vec3d &cPos = mModel->mOpenedC3dFile->GetMarkerPos(f, i);
		if ( !(cPos[0] == 0 && cPos[1] == 0 && cPos[2] == 0) ) {	//Skip (0,0,0)
			constraintIDs.push_back(i);
		}
	}
}

/**
* Builds a Jacobian Matrix J[]
*/
Matd Solver::computeJ(Markers handles, int cID){
	Matd J;
	J.SetSize(mModel->GetDofCount(), 4);
	J.MakeZero();

	//Get the current Node
	Marker *h = handles[cID];
	Vec4d tempPos(h->mOffset, 1.0);
	TransformNode *node = mModel->mLimbs[h->mNodeIndex];

	// Utility identity matrix
	Mat4d identity;
	identity.MakeDiag(1.0);

	// Recursively create Jacobian
	computeJ(handles, cID, J, node, tempPos, identity);

	return J;
}

/*
* Jacobian recursive worker method
*/
void computeJ(std::vector<Marker*> handles, int cID, Matd &J, TransformNode *node, Vec4d &pos, Mat4d &identity) {

	// Get initial transform data
	std::vector<Transform *> transforms = node->mTransforms;

	// loop over all transforms in node, figuring out
	// appropriate pre-matrix for each DOF transform
	Mat4d leftMat;
	preMatrix.MakeDiag(1.0);	// start off as identity
	for (int i = 0; i < transforms.size(); i++)
	{
		Transform * transform = transforms[i];
		// only need to calculate if DOF transform
		if (transform->IsDof())
		{
			// loop through each dof
			for (int j = 0; j < transform->GetDofCount(); j++)
			{
				// get dof
				Dof * dof = transform->GetDof(j);
				int dofId = dof->mId;

				// place in appropriate spot in map
				mPreMatrices[dofId] = preMatrix;
			}
		}
		// now that we've gone through this transform, 
		// update the preMatrix to have current transform for later transforms
		preMatrix = preMatrix * transform->GetTransform();
	}
}




