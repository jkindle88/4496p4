#ifndef __COMMAND_H__
#define __COMMAND_H__

#include <fstream>
#include "vl/VLd.h"
#include <vector>
#include "Transform.h"
#include "TransformNode.h"
#include "RotateEuler.h"
#include "Translate.h"

typedef void (*Command)(void*);

void LoadModel(void*);
void LoadC3d(void*);
void Exit(void*);
void Solution(void*);
void computeJ(int dofIndex, int markIndex, TransformNode *t, Mat4d parentT, Matd &J);
#endif

