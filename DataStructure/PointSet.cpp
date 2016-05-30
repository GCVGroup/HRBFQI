// PointSet.cpp: PointSet クラスのインプリメンテーション
//
//////////////////////////////////////////////////////////////////////
#include <afxwin.h>         // MFC のコアおよび標準コンポーネント
#include <afxext.h>         // MFC の拡張部分
#include <afxdisp.h>        // MFC のオートメーション クラス
#include <afxdtctl.h>		// MFC の Internet Explorer 4 コモン コントロール サポート
#ifndef _AFX_NO_AFXCMN_SUPPORT
#include <afxcmn.h>			// MFC の Windows コモン コントロール サポート
#endif // _AFX_NO_AFXCMN_SUPPORT

#include "PointSet.h"
#include "math.h"
#include "../numericalC/jacobi.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

PointSet::PointSet()
{
	point = NULL;
	normal = NULL;
	color = NULL;
	value = NULL;
	quad = NULL;
}
PointSet::PointSet(int N, float (*_point)[3], float (*_normal)[3])
{
	setPointSize(N);
	for(int i=0; i<N; i++){
	  setPoint(i, _point[i][0], _point[i][1], _point[i][2]);
	  setNormal(i, _normal[i][0], _normal[i][1], _normal[i][2]);
	}
}
PointSet::~PointSet()
{
	if(point != NULL)
		delete[] point;
	if(normal != NULL)
		delete[] normal;
	if(color != NULL)
		delete[] color;
	if(value != NULL)
		delete[] value;
	if(quad != NULL){
		delete[] quad;
		delete[] tangent1;
		delete[] tangent2;

	}
}
void PointSet::rescale(float dia)
{
	float max[3], min[3];
	getBound(min, max, 1.0f);
	float sx = max[0] - min[0];
	float sy = max[1] - min[1];
	float sz = max[2] - min[2];
	float diaO = (float)sqrt(sx*sx + sy*sy + sz*sz);
	float scale = dia/diaO;
	for(int i=0; i<point_N; i++)
	{
		float* p = point[i];
		p[0] *= scale;
		p[1] *= scale;
		p[2] *= scale;
	}
}
void PointSet::scale(float ori[3], float scale)
{
	int i;
	for(i=0; i<point_N; i++){
		float* p = point[i];
		p[0] = (p[0]-ori[0])*scale;
		p[1] = (p[1]-ori[1])*scale;
		p[2] = (p[2]-ori[2])*scale;
	}
}
float PointSet::fitIntoBox(float ct[3], float boxSize)
{
	float xMin, xMax, yMin, yMax, zMin, zMax;
	getBound(xMin, xMax, yMin, yMax, zMin, zMax);
	//  calculate scale to fit into a 2*2*2 box
//	printf("Original bounding box: (%f, %f, %f)-(%f, %f, %f)\n", xMin, yMin, zMin, xMax, yMax, zMax);
//	return 0;

	float sX = 	2.0/(xMax - xMin);
	float sY = 	2.0/(yMax - yMin);
	float sZ = 	2.0/(zMax - zMin);
	float scale = sX;
	if (sY < scale)
	  scale = sY;
	if (sZ < scale)
	  scale = sZ;
	scale *= boxSize;

	float cX = (xMax + xMin) / 2.0f;
	float cY = (yMax + yMin) / 2.0f;
	float cZ = (zMax + zMin) / 2.0f;
	
	ct[0] = cX;	ct[1] = cY;	ct[2] = cZ;
	int numVert = point_N;
	for (int i = 0; i < numVert; i++){
	//  pointer to vertex information
	float *p = point[i];
	p[0] = (p[0]-cX)*scale;
	p[1] = (p[1]-cY)*scale;
	p[2] = (p[2]-cZ)*scale;
	}

	return scale;
}

void PointSet::swapIndex(int i, int j)
{
	float tmp = point[i][0];
	point[i][0] = point[j][0];
	point[j][0] = tmp;

	tmp = point[i][1];
	point[i][1] = point[j][1];
	point[j][1] = tmp;

	tmp = point[i][2];
	point[i][2] = point[j][2];
	point[j][2] = tmp;

	if(normal != NULL){
		tmp = normal[i][0];
		normal[i][0] = normal[j][0];
		normal[j][0] = tmp;

		tmp = normal[i][1];
		normal[i][1] = normal[j][1];
		normal[j][1] = tmp;

		tmp = normal[i][2];
		normal[i][2] = normal[j][2];
		normal[j][2] = tmp;
	}

	if(value != NULL){
		tmp = value[i];
		value[i] = value[j];
		value[j] = tmp;
	}
}

void PointSet::setPointSize(int n)
{
	this->point_N = n;
	point = new float[n][3];
	normal = new float[n][3];
}

void PointSet::setPoint(int i, float x, float y, float z)
{
	point[i][0] = x;
	point[i][1] = y;
	point[i][2] = z;
}

void PointSet::setNormal(int i, float x, float y, float z)
{
    if(normal == NULL)
      normal = new float[point_N][3];
	normal[i][0] = x;
	normal[i][1] = y;
	normal[i][2] = z;
}
void PointSet::setColor(int i, float r, float g, float b)
{
	if(color == NULL)
		color = new float[point_N][3];
	color[i][0] = r;
	color[i][1] = g;
	color[i][2] = b;
}

void PointSet::setValue(int i, float v)
{
	if(value == NULL)
		value = new float[point_N];
	value[i] = v;
}

void PointSet::getBound(float &xmin, float &xmax, 
						 float &ymin, float &ymax, 
						 float &zmin, float &zmax)
{
	xmin = point[0][0];
	xmax = point[0][0];
	ymin = point[0][1];
	ymax = point[0][1];
	zmin = point[0][2];
	zmax = point[0][2];

	for(int i=1; i<point_N; i++){
		if(point[i][0] > xmax)
			xmax = point[i][0];
		if(point[i][0] < xmin)
			xmin = point[i][0];

		if(point[i][1] > ymax)
			ymax = point[i][1];
		if(point[i][1] < ymin)
			ymin = point[i][1];

		if(point[i][2] > zmax)
			zmax = point[i][2];
		if(point[i][2] < zmin)
			zmin = point[i][2];
	}
}
void PointSet::getBound(float min[3], float max[3], float rate)
{
	for(int i=0; i<3; i++)
	{
		max[i] = min[i] = point[0][i];
		for(int j=1; j<point_N; j++)
		{
			if(min[i] > point[j][i])
			min[i] = point[j][i];
			if(max[i] < point[j][i])
			max[i] = point[j][i];
		}
		float s = max[i]-min[i];
		float m = 0.5f*(max[i]+min[i]);
		max[i] = m + 0.5f*rate*s;
		min[i] = m - 0.5f*rate*s;
	}
}

void PointSet::getBound(float min[3], float max[3], int start, int end)
{
	for(int i=0; i<3; i++)
	{
		max[i] = min[i] = point[start][i];
		for(int j=start+1; j<end; j++)
		{
			if(min[i] > point[j][i])
			min[i] = point[j][i];
			if(max[i] < point[j][i])
			max[i] = point[j][i];
		}
		float s = max[i]-min[i];
		float m = 0.5f*(max[i]+min[i]);
		max[i] = m + 0.5f*s;
		min[i] = m - 0.5f*s;
	}
}

//Index "end" is not taken int account
void PointSet::centroid(float c[3], int start, int end)
{
	c[0] = c[1] = c[2] = 0;
	for(int i=start; i<end; i++)
	{
		c[0] += point[i][0];
		c[1] += point[i][1];
		c[2] += point[i][2];
	}
	c[0] /= (end-start);
	c[1] /= (end-start);
	c[2] /= (end-start);
}

//Index "end" is not taken int account
void PointSet::averagedNormal(float n[3], int start, int end)
{
	n[0] = n[1] = n[2] = 0;
	for(int i=start; i<end; i++)
	{
		n[0] += normal[i][0];
		n[1] += normal[i][1];
		n[2] += normal[i][2];
	}
	double len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
	if((float)len != 0)
	{
		n[0] = (float)(n[0]/len);
		n[1] = (float)(n[1]/len);
		n[2] = (float)(n[2]/len);
	}
}

void PointSet::computeNormalWithCV(float n[3], int* list, int N, PointSet *pts)
{
    int i;
	PointSet *tempPs = this;
	if(pts != 0)
		tempPs = pts;

    float (*tempPoint)[3] = tempPs->point;
    
    float cx = 0;
    float cy = 0;
    float cz = 0;
    
    for(i=0; i<N; i++){
      float* p = tempPoint[list[i]];
      cx += p[0];
      cy += p[1];
      cz += p[2];
    }
    float Ni = 1.0f/N;
    cx *= Ni;
    cy *= Ni;
    cz *= Ni;
    
    //data for Jacobi method
    float **A = new float*[4];
    float **v = new float*[4];
    float w[4];
    int nrot;
    for(i=1; i<4; i++){
      A[i] = new float[4];
      A[i][1] = A[i][2] = A[i][3] = 0;
      v[i] = new float[4];
    }
    
    //CV matrix
    for(i=0; i<N; i++){
      float* p = tempPoint[list[i]];
      
      float vx = p[0] - cx;
      float vy = p[1] - cy;
      float vz = p[2] - cz;
      
      A[1][1] += vx*vx;
      A[1][2] += vx*vy;
      A[1][3] += vx*vz;
      
      A[2][2] += vy*vy;
      A[2][3] += vy*vz;
      
      A[3][3] += vz*vz;
    }
    A[2][1] = A[1][2];
    A[3][1] = A[1][3];
    A[3][2] = A[3][2];
    
    Jacobi jacobiSolver;
	jacobiSolver.jacobi(A, 3, w, v, &nrot);
    
    int mini;
    if(fabs(w[1]) < fabs(w[2]))
      mini = 1;
    else
      mini = 2;
    if(fabs(w[mini]) > fabs(w[3]))
      mini = 3;
    
    n[0] = v[1][mini];
    n[1] = v[2][mini];
    n[2] = v[3][mini];
    
    for(i=1; i<4; i++){
      delete[] A[i];
      delete[] v[i];
    }
    delete[] A;
    delete[] v;	
}
