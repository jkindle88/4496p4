/**
* Solver.h
* @Author: Joe Kindle
* Worker class for solving inverse kinematics.
* Will also be used to store the relevant values.
*/

//Includes
#include <map>
#include "Model.h"

class Solver
{
public:
	Solver(double epsilon, double step, int maxIterations, int maxFrames, Model *model);
	~Solver();

	//Vars
	double mEps;
	double mStep;
	int mMaxIters;
	int mMaxFrames;
	Model *mModel;
	vector<int> constraintIDs;
	vector<Vec3d> constraintVals;
	vector<double> constraintLengths;

	// Used by Jacobian computation
	// Attempting to speed up lookup a little bit...
	// Newton-Raphson method
	// Used to map an ID to a Matrix
	std::map<int, Mat4d> left;
	std::map<int, Mat4d> right;

	//Main Loop
	void solve();
	//Compute Constraint Matrix
	void buildConstraintMat(int frameNumber);
	//Jacobian fun
	Matd computeJ(std::vector<Marker*> handles, int cID);
	void computeJ(std::vector<Marker*> handles, int cID, Matd &J, TransformNode *node, Vec4d &pos, Mat4d &identity);
};

//Type definition for the current marker handles
typedef std::vector<Marker*> Markers;
