/**
* Solver.h
* @Author: Joe Kindle
* Worker class for solving inverse kinematics.
* Will also be used to store the relevant values.
*/

//Includes
#include "Model.h"

class Solver
{
public:
	Solver(double epsilon, double step, int maxIterations, int maxFrames, Model *model);
	~Solver();

	//Main Loop
	void solve();

	//Utility functions to calculate constraints and gradients
	void computeConstraints(int frameNumber,bool create);
	Vecd computeGradient(int frameNumber);

	//Vars
	double mEps;
	double mStep;
	int mMaxIters;
	int mMaxFrames;
	Model *mModel;

};

//Type definition for the current marker handles
typedef std::vector<Marker*> Markers;