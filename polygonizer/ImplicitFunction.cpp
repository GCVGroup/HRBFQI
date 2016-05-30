// ImplicitFunction.cpp: 
//
//////////////////////////////////////////////////////////////////////

#include "ImplicitFunction.h"

//////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////

ImplicitFunction::ImplicitFunction()
{

}

ImplicitFunction::~ImplicitFunction()
{

}

float ImplicitFunction::value(float x, float y, float z)
{
	return 0;
}

void ImplicitFunction::gradient(float g[], float x, float y, float z)
{

}

float ImplicitFunction::value(float x, float y, float z, bool &isValid)
{
	return 0;
}

void ImplicitFunction::asignValueToVoxels(float ***&v, bool ***&isIn, float o[], float space[], int dim[])
{

}
