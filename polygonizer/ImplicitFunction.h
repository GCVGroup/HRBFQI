// ImplicitFunction.h: ImplicitFunction 
//////////////////////////////////////////////////////////////////////

//#if !defined(AFX_IMPLICITFUNCTION_H__64F070AD_4EFC_4106_93D7_3283C4A08551__INCLUDED_)
//#define AFX_IMPLICITFUNCTION_H__64F070AD_4EFC_4106_93D7_3283C4A08551__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class ImplicitFunction  
{
public:
	virtual void asignValueToVoxels(float ***&v, bool ***&isIn, float o[], float space[], int dim[]);
	virtual float value(float x, float y, float z, bool &isValid);
	virtual void gradient(float g[3], float x, float y, float z);
	virtual float value(float x, float y, float z);
	ImplicitFunction();
	virtual ~ImplicitFunction();

};

//#endif // !defined(AFX_IMPLICITFUNCTION_H__64F070AD_4EFC_4106_93D7_3283C4A08551__INCLUDED_)
