#ifndef __COMMAND_H__
#include "Command.h"
#endif //__COMMAND_H__

#ifndef __C3DFILEINFO_H__
#include "C3dFileInfo.h"
#endif	//__C3DFILEINFO_H__

#ifndef __ARTICULATEDBODY_H__
#include "ArticulatedBody.h"
#endif	//__ARTICULATEDBODY_H__

#ifndef RealTimeIKui_h
#include "RealTimeIKui.h"
#endif //RealTimeIKui_h

#ifndef __PHYLTERFLGLWINDOW_H__
#include "PhylterGLWindow.h"
#endif	//__PHYLTERFLGLWINDOW_H__

#ifndef	__PHOWARDDATA_H__
#include "PhowardData.h"
#endif

#ifndef __TRANSFORM_H__
#include "Transform.h"
#endif	//__TRANSFORM_H__

#include "Solver.h"



int readSkelFile( FILE* file, ArticulatedBody* skel );

extern RealTimeIKUI *UI;

void LoadModel(void *v)
{
  char *params = (char*)v;
  if(!params){
    params = (char*)fl_file_chooser("Open File?", "{*.skel}", "../4496p4/src/skels" );
  }

  if(!params)
    return;

  FILE *file = fopen(params, "r");
    
  if(file == NULL){
    cout << "Skel file does not exist" << endl;
    return;
  }

  ArticulatedBody *mod = new ArticulatedBody();
  UI->mData->mModels.push_back(mod);
  UI->mData->mSelectedModel = mod;

  readSkelFile(file, mod);
  UI->CreateDofSliderWindow();

  mod->InitModel();
  UI->mGLWindow->mShowModel = true;
  UI->mShowModel_but->value(1);
  UI->mGLWindow->refresh();
  
  cout << "number of dofs in model: " << UI->mData->mModels[0]->GetDofCount() << endl;
}

void Exit(void *v)
{
  exit(0);
}

void LoadC3d(void *v)
{
  if(!UI->mData->mSelectedModel){
    cout << "Load skeleton first";
    return;
  }
  char *params = (char*)v;
  if(!params){
    params = fl_file_chooser("Open File?", "{*.c3d}", "../4496p4/src/mocap" );
  }

  if(!params)
    return;
  
  char *c3dFilename = new char[80];
  
  // load single c3d file
 
  C3dFileInfo *openFile = new C3dFileInfo(params);
  openFile->LoadFile();
  UI->mData->mSelectedModel->mOpenedC3dFile = openFile;
  cout << "number of frames in c3d: " << openFile->GetFrameCount() << endl;

  UI->InitControlPanel();
  UI->mGLWindow->mShowConstraints = true;
  UI->mShowConstr_but->value(1);
}

void Solution(void *v)
{
    cout << "TODO: Solve inverse kinematics problem" << endl;
    bool test = UI->mData->mSelectedModel->mLimbs[0]->mTransforms[0]->IsDof();
	
	double eps = DBL_EPSILON; // minimize the step size
	double step = DBL_EPSILON;
	int maxIterations = 3;
	int maxFrames = 100; //play with this

	//Store the model 
	Model *m = UI->mData->mSelectedModel;
	//Store the list of marker handles on the model
	Markers &handles = m->mHandleList;

    int numFrames = m->mOpenedC3dFile->GetFrameCount();
    int numDof = m->GetDofCount();
    int numMarkers = m->GetHandleCount();
    int numLimbs = m->GetNodeCount();

    Matd frames = Matd(numFrames, numDof);

	for(int f = 0; f < numFrames && f < maxFrames; f++)
	{
		bool done = false;

		while(!done)
		{
			done = false;
			//Build the relevant matricies
			Vecd q = Vecd(numDof);
			Vecd C = Vecd(numMarkers*3);
			Matd J = Matd(numMarkers*3, numDof, vl_0);

			//Store the list of marker handles on the model
			Markers &handles = m->mHandleList;

			m->mDofList.GetDofs(&frames[f]);

			//Build C() constraint matrix
			for(int i = 0; i < numMarkers*3; i+=3)
			{
					Vec3d markerPos = m->mOpenedC3dFile->GetMarkerPos(f, i/3);

					Vec3d c = (markerPos - m->mHandleList[i/3]->mGlobalPos);
					C[i*3] = c[0];
					C[i*3+1] = c[1];
					C[i*3+2] = c[2];

					//Halting condition
					if(abs(c[0]) > eps || abs(c[1]) > eps || abs(c[2]) > eps)
					{
							done = true;
					}
			}

			//Build J() Jacobian
			for(int markI = 0; markI < numMarkers; markI++)
			{
					for(int dofI = 0; dofI < numDof; dofI++)
					{
							int limbI = numLimbs-1;
							SetJacobianEntry(dofI, markI, m->mLimbs[limbI], Mat4d(vl_one), J); 
					}
			}

			//compute J*Jtranspose
			Matd JJtrans = J*trans(J);
			//take the inverse of J*Jtranspose
			Matd JJtInv = inv(JJtrans);
			//take the Moore-Penrose PI
			Matd MPPI_J = trans(J)*JJtInv;

			//Add current step
			frames[f] = frames[f] + (MPPI_J*C);

			//Store results
			m->SetDofs(frames, f);
		}
	}

	//Create and run the solver
	//Solver s(eps, step, maxIterations, maxFrames, UI->mData->mSelectedModel);

	//s.solve();

}

void SetJacobianEntry(int dofIndex, int markIndex, TransformNode *t, Mat4d parentTransform, Matd &Jacobian)
{
        /*start with the product of all parent transforms*/
        Mat4d currentTransform = parentTransform;

        /*Apply transforms associated with this node
         * If this node contains the Dof, apply the derivative of that tranform*/
        for(int transformIndex = 0; transformIndex < t->mTransforms.size(); transformIndex++)
        {       
                if(t->mTransforms[transformIndex]->IsDof())
                {
                        bool containsTHEDof = false;
                        for(int dof = 0; dof < t->mTransforms[transformIndex]->GetDofCount(); dof++)
                        {
                                if(t->mTransforms[transformIndex]->GetDof(dof)->mId == dofIndex)
                                {
                                        Mat4d derivativeTransform = t->mTransforms[transformIndex]->GetDerivative(0,dof);
                                        currentTransform = currentTransform * derivativeTransform;
                                        containsTHEDof = true;
                                }
                        }
                        if(!containsTHEDof)
                        {
                                currentTransform = currentTransform * t->mTransforms[transformIndex]->GetTransform(0);
                        }
                }
                else
                {
                        currentTransform = currentTransform * t->mTransforms[transformIndex]->GetTransform(0);
                }
        }

        /*If this node holds the marker we are evaluating at,
         * Multiply by that offset and throw the result into the Jacobian*/
        for(int i = 0; i < t->mHandles.size(); i++)
        {
                if(t->mHandles[i]->mNodeIndex == markIndex)
                {
                        Vec4d hk = Vec4d(t->mHandles[i]->mOffset, 1);
                        Vecd JacobianEntry = currentTransform * hk;

                        double x[3];
                        x[0] = JacobianEntry[0];
                        x[1] = JacobianEntry[1];
                        x[2] = JacobianEntry[2];

                        Jacobian[(markIndex*3)+0][dofIndex] = JacobianEntry[0];
                        Jacobian[(markIndex*3)+1][dofIndex] = JacobianEntry[1];
                        Jacobian[(markIndex*3)+2][dofIndex] = JacobianEntry[2];

                        return;
                }
        }

        /*Otherwise, continue on to the children*/
        for(int childIndex = 0; childIndex < root->mChildren.size(); childIndex++)
        {
                SetJacobianEntry(handleId, dofId, ((TransformNode*)root->mChildren[childIndex]), currentTransform, Jacobian);
        }
}

