// PointSet.h: PointSet 
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_POINTSET_H__E3B75762_65C7_489D_A051_4FE7277A5630__INCLUDED_)
#define AFX_POINTSET_H__E3B75762_65C7_489D_A051_4FE7277A5630__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class PointSet  
{
public:
	PointSet();
	PointSet(int N, float (*_point)[3], float (*_normal)[3]);
	virtual ~PointSet();

public:
	void computeNormalWithCV(float n[3], int* list, int N, PointSet *pts = 0);
	void swapIndex(int i, int j);
	void rescale(float dia);
	float fitIntoBox(float ct[3], float boxSize = 10.0f);
	void scale(float ori[3], float scale);
	void getBound(float min[3], float max[3], float rate);
	void getBound(float min[3], float max[3], int start, int end);
	void centroid(float c[3], int start, int end);
	void averagedNormal(float n[3], int start, int end);
	void getBound(float &xmin, float &xmax, float &ymin, float &ymax, float &zmin, float &zmax);
	void setValue(int i, float v);
	void setNormal(int i, float x, float y, float z);
	void setPoint(int i, float x, float y, float z);
	void setColor(int i, float r, float g, float b);
	void setPointSize(int n);
	float (*point)[3];
	float (*normal)[3];
	float (*color)[3];

	float *value;

	float (*quad)[3];
	float (*tangent1)[3];
	float (*tangent2)[3];
	int point_N;
};

#endif // !defined(AFX_POINTSET_H__E3B75762_65C7_489D_A051_4FE7277A5630__INCLUDED_)
