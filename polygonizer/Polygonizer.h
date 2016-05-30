// Polygonizer.h: Polygonizer 
//
//////////////////////////////////////////////////////////////////////

#include "../DataStructure/PointSet.h"
#include "../DataStructure/PolygonalMesh.h"
#include "../numericalC/SVD.h"
#include "ImplicitFunction.h"

class Polygonizer  
{
public:
	float spaceX, spaceY, spaceZ;
	int dimX, dimY, dimZ;
	float originX, originY, originZ;
	ImplicitFunction* func;

	int thread_num;
public:
	void searchZero(float p[], float start[], float end[], float f1, float f2, float e);
	void cutQuad(PolygonalMesh* mesh);
	void bisection(float p[3], float start[3], float end[3], float e);
	PolygonalMesh* dualContouring(float epsilon, float tau);
	float smoothGridPointValue(int i, int j, int k, float ***fValue, bool ***pValid);
	bool isBoundaryGridPoint(int i, int j, int k, float ***fValue, bool ***pValid);

	inline void weightD(float g[], double d2, double R2, float vx, float vy, float vz);
	inline float value(float vx, float vy, float vz, float nx, float ny, float nz, double R2, bool &isValid);

	void setSpace(float x, float y, float z);
	void setOrigin(float x, float y, float z);
	void setDim(int dimX, int dimY, int dimZ);
	Polygonizer();
	virtual ~Polygonizer();


	/*
	static inline BOOL INVERSE(double B[10], double A[10]){
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
	*/

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

	static inline void MAT_SUM(double B[10], double A[10]){
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
