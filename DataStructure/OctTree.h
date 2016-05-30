// OctTree.h: OctTree 
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_OCTTREE_H__B2DDC163_33E0_4996_8D78_7153D535816E__INCLUDED_)
#define AFX_OCTTREE_H__B2DDC163_33E0_4996_8D78_7153D535816E__INCLUDED_

//#if _MSC_VER > 1000
//#pragma once
//#endif // _MSC_VER > 1000

#include "PointSet.h"
#include "math.h"

#define MAX 8
#define EPSILON 0.000001

class OctTree  
{
	class Cell{
		float ox;
		float oy;
		float oz;
		int level;
		Cell* child[8];
		bool isLeaf;
		int* index;
		int index_N;
		float p[3], n[3];

		public:
			void addDoubledPointIntoList(PointSet *ps, int &size, float (*&p_list)[3], float *&value, float off);
			void cleanChild();
			void splitChild(int M, PointSet* ps, float sizeX, float sizeY, float sizeZ);
			void computeCellPoint(PointSet* ps, double *&Q);
			void addPointIntoList(PointSet* ps, int &size, float **p_list, float **n_list, float x, float y, float z, float r, int max);
			void normalizeNormals();
			void computeCenter(PointSet* ps);
			void addLeafSize(float sizeX, float sizeY, float sizeZ, int &count, float &sum);
			void addIntoTable(PointSet *ps, int* table, int &tableSize, float x, float y, float z, float r, float sizeX, float sizeY, float sizeZ);
			void createChild(float sizeX, float sizeY, float sizeZ);
			void addPoint(int i, PointSet* ps, float sizeX, float sizeY, float sizeZ);
			inline float getSizeZ(float sizeZ);
			inline float getSizeY(float sizeY);
			inline float getSizeX(float sizeX);
			Cell(int level, float ox, float oy, float oz);
			virtual ~Cell();
	};

public:
	OctTree();
	virtual ~OctTree();

private:

public:
	void setDoubledPoints(int &size, float (*&p_list)[3], float *&value);
	void getPointList(int &size, float **&p_list, float **&n_list, float x, float y, float z, float r, int max);
	void computeCenter();
	float getAverageLeafSize();
	int* getIndexTable(int &tableSize, float x, float y, float z, float r);
	bool getIndexTable(int *table, int &tableSize, float x, float y, float z, float r);
	void setPointSet(PointSet *ps);
	PointSet* ps;

	int max;

	float originX;
	float originY;
	float originZ;

	float sizeX;
	float sizeY;
	float sizeZ;

	Cell* root;

private:
	int* table;
	float **point_list;
	float **normal_list;


public:
	void splitChild(int M);
	void computeCellPoint();
	static inline bool INVERSE(double B[10], double A[10]){
		double d = DET(A);
		if(fabs(d) < EPSILON)
			return false;
		B[0] = (A[3]*A[5] - A[4]*A[4])/d;
		B[1] = (A[2]*A[4] - A[1]*A[5])/d;
		B[2] = (A[1]*A[4] - A[2]*A[3])/d;
		B[3] = (A[0]*A[5] - A[2]*A[2])/d;
		B[4] = (A[1]*A[2] - A[0]*A[4])/d;
		B[5] = (A[0]*A[3] - A[1]*A[1])/d;
		return true;
	}

	static inline double DET(double A[10]){
		return A[0]*A[3]*A[5] + 2.0*A[1]*A[4]*A[2] 
			-A[2]*A[2]*A[3] - A[1]*A[1]*A[5] - A[4]*A[4]*A[0];
	}

	static inline void MATRIX(double A[10], double n[3], double d){
		A[0] = n[0]*n[0];
		A[1] = n[0]*n[1];
		A[2] = n[0]*n[2];
		A[3] = n[1]*n[1];
		A[4] = n[1]*n[2];
		A[5] = n[2]*n[2];
		A[6] = d*n[0];
		A[7] = d*n[1];
		A[8] = d*n[2];
		A[9] = d*d;
	}

	static inline void MAT_TIMES(double A[10], double k){
		A[0] *= k;
		A[1] *= k;
		A[2] *= k;
		A[3] *= k;
		A[4] *= k;
		A[5] *= k;
		A[6] *= k;
		A[7] *= k;
		A[8] *= k;
		A[9] *= k;
	}

	static inline void MAT_SUM(double B[6], double A[6]){
		B[0] += A[0];
		B[1] += A[1];
		B[2] += A[2];
		B[3] += A[3];
		B[4] += A[4];
		B[5] += A[5];
		B[6] += A[6];
		B[7] += A[7];
		B[8] += A[8];
		B[9] += A[9];
	}

	static inline void MAT_BY_VEC(double v[3], double A[10], double b[3]){
		v[0] = A[0]*b[0] + A[1]*b[1] + A[2]*b[2];
		v[1] = A[1]*b[0] + A[3]*b[1] + A[4]*b[2];
		v[2] = A[2]*b[0] + A[4]*b[1] + A[5]*b[2];
	}

	static inline void MAT_BY_VEC(float v[3], double A[10], float b[3]){
		v[0] = (float)(A[0]*b[0] + A[1]*b[1] + A[2]*b[2]);
		v[1] = (float)(A[1]*b[0] + A[3]*b[1] + A[4]*b[2]);
		v[2] = (float)(A[2]*b[0] + A[4]*b[1] + A[5]*b[2]);
	}

	static inline void MAT_INIT(double A[10]){
		A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = A[9] = 0;
	}

	static inline void MAT_PLUS(double C[10], double A[10], double B[10]){
		C[0] = A[0] + B[0];
		C[1] = A[1] + B[1];
		C[2] = A[2] + B[2];
		C[3] = A[3] + B[3];
		C[4] = A[4] + B[4];
		C[5] = A[5] + B[5];
		C[6] = A[6] + B[6];
		C[7] = A[7] + B[7];
		C[8] = A[8] + B[8];
		C[9] = A[9] + B[9];
	}

	static inline void MAT_COPY(double A[10], double B[10]){
		A[0] = B[0];
		A[1] = B[1];
		A[2] = B[2];
		A[3] = B[3];
		A[4] = B[4];
		A[5] = B[5];
		A[6] = B[6];
		A[7] = B[7];
		A[8] = B[8];
		A[9] = B[9];
	}

	static inline double Q_ERR(double A[10], float v[3]){
		return v[0]*(A[0]*v[0] + A[1]*v[1] + A[2]*v[2]) +
			   v[1]*(A[1]*v[0] + A[3]*v[1] + A[4]*v[2]) +
			   v[2]*(A[2]*v[0] + A[4]*v[1] + A[5]*v[2]) +
			   2.0*(A[6]*v[0] + A[7]*v[1] + A[8]*v[2]) + 
			   A[9];
	}

	static inline double Q_ERR(double A[10], double v[3]){
		return v[0]*(A[0]*v[0] + A[1]*v[1] + A[2]*v[2]) +
			   v[1]*(A[1]*v[0] + A[3]*v[1] + A[4]*v[2]) +
			   v[2]*(A[2]*v[0] + A[4]*v[1] + A[5]*v[2]) +
			   2.0*(A[6]*v[0] + A[7]*v[1] + A[8]*v[2]) + 
			   A[9];
	}
};

#endif // !defined(AFX_OCTTREE_H__B2DDC163_33E0_4996_8D78_7153D535816E__INCLUDED_)
