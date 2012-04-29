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
	Solver::leftMap;
	Solver::rightMap;
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
	const int sIFreq = 20;
	const double sIF = 4.0;
	const double sDF = 2.0;
	const int eIFreq = 2;
	const double eIF = 2.0;
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

				//Calculate current and add to gradient matrix
				Vecd currentVal = J * cVec;
				gradient = gradient + currentVal;

			}

			// Finally, double it
			gradient = 2.0 * gradient;

			// Get previous DOFs
			Vecd prevDofs;
			prevDofs.SetSize(mModel->GetDofCount());
			mModel->mDofList.GetDofs(&prevDofs);

			// Update DOFs
			Vecd nextDofs = prevDofs - s * gradient;
			mModel->SetDofs(nextDofs);

			
			// Recalculate constraints and objective function
			double newO = 0.0;

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
				newO += cLength; 
			}


			// Attempting to make sure we don't go infinite...
			if (newO < o)
			{		
				if (newO > mEps && iterCount > 0 && iterCount % 20 == 0)
				{
					s *= 3.0;
				}
				o = newO;	
			}
			else
			{
				stepCount++;
				// Decrease step size
				s = s / 2.0;
				
				// based on how many times we've had to decrease the step size, choose different options
				if (stepCount < mMaxIters)
				{
					// Reset to the previous dofs and try again
					mModel->SetDofs(prevDofs);
				}
				else if (stepCount < eIFreq*mMaxIters)
				{
					// Try increasing epsilon
					e *= 3.0;

					// Reset to the previous dofs and try again
					mModel->SetDofs(prevDofs);
				}
				else
				{
					// Otherwise, halt or we will never finish
					o = 0.0;
				}
			}
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
	computeJ(J, node, tempPos, identity);

	return J;
}

/*
* Jacobian recursive worker method
*/
void Solver::computeJ(Matd &J, TransformNode *node, Vec4d &pos, Mat4d &trans) {

	// Get initial transform data
	std::vector<Transform *> transforms = node->mTransforms;

	// LEFT MATRIX
	Mat4d leftMat;
	leftMat.MakeDiag(1.0);	// identity
	for (int i = 0; i < transforms.size(); i++)
	{
		Transform *tr = transforms[i];
		// Only need to calculate if DOF
		if (tr->IsDof())
		{
			for (int j = 0; j < tr->GetDofCount(); j++)
			{
				// Get DOF
				Dof *dof = tr->GetDof(j);
				int dofId = dof->mId;

				// Add to our map
				leftMap[dofId] = leftMat;
			}
		}
		// Now update leftMat
		leftMat = leftMat * tr->GetTransform();
	}

	// RIGHT MATRIX
	Mat4d rightMat;
	rightMat.MakeDiag(1.0);	// identity
	for (int i = transforms.size()-1; i >= 0; i--)
	{
		Transform *tr = transforms[i];
		// Only need to calculate if DOF
		if (tr->IsDof())
		{
			for (int j = tr->GetDofCount()- 1; j >= 0 ; j--)
			{
				// Get DOF
				Dof *dof = tr->GetDof(j);
				int dofId = dof->mId;

				// Add to our map
				rightMap[dofId] = rightMat;
			}
		}
		// Now update leftMat
		rightMat = tr->GetTransform() * rightMat;
	}

	// JACOBIAN
	Mat4d pTrans = node->mParentTransform;

	// New identity matrix for later use
	Mat4d newTransform;
	newTransform.MakeDiag(1.0);

	for (int i = 0; i < transforms.size(); i++)
	{		
		// Check if DOF, if so compute derivative
		Transform *tr = transforms[i];
		if (tr->IsDof())
		{
			for (int j = 0; j < tr->GetDofCount(); j++)
			{
				Mat4d deriv = tr->GetDeriv(j);

				// Get row
				Dof * dof = tr->GetDof(j);
				int dofId = dof->mId;

				Mat4d leftMat = leftMap[dofId];
				Mat4d rightMat = rightMap[dofId];
				Vec4d value = pTrans * leftMat * deriv * rightMat * trans * pos;

				J[dofId] = value;
			}
		}

		newTransform = newTransform * tr->GetTransform();
	}
	// Calculate J[] values for parent if necessary
	TransformNode *parent = node->mParentNode;
	// If the parent is non-null and not the current node
	if (parent != NULL && parent != node)
	{
		newTransform = newTransform * trans;
		computeJ(J, parent, pos, newTransform);
	}
}




