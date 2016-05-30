#pragma once
#include "../DataStructure/PointSet.h"
#include "../DataStructure/OctTree.h"
#include "math.h"
#include "../DataStructure/PolygonalMesh.h"
#include "../numericalC/SVD.h"
#include "../polygonizer/ImplicitFunction.h"

class HRBF :
	public ImplicitFunction
{
public:
	PointSet* ps;
	OctTree* tree;
	float support;

	//Shift amount
	double *sol;						

	//	------------------------------------
	int maxNeighborNumber;

private:
	//T2 = support^2
	double T2;

	int thread_num;

public:
	HRBF(void);
	~HRBF(void);

	void differenceFuncGradientAndNormalPt(double &maxAng, double &aveAng);

	float value(float x, float y, float z);
	float value(float x, float y, float z, bool &isValid);
	void gradient(float g[], float x, float y, float z);
	float valueAndGradient(float g[], float x, float y, float z);

	void setSupport(float T);
	int getMaximalNeighborsInSupport(float T);

	//Octree oprtations
	void setPointSet(PointSet* ps);
	float getAveragedLeafSize();

	//Wendland's CSRBF
	inline void weightD(float *w, float g[], double d2, float vx, float vy, float vz);
	inline void weightD(float g[], double d2, float vx, float vy, float vz);
	inline void weightH(float h[], double d2, float vx, float vy, float vz);
	inline void weightDH(float g[], float h[], double d2, float vx, float vy, float vz);
	inline void weightDH(float *w, float g[], float h[], double d2, float vx, float vy, float vz);

	inline double weight(double d2);

	//Visualizations
	PolygonalMesh* generateCrossSection(float o[3], float t1[3], float t2[3], int n, int m);

	//Fitting
	void fit(float support, float nsmooth);

	//	Implementation of the paper "Hermite Radial Basis Functions Implicits"
	//	single level
	void computeSolution(float nsmooth);

	//	quasi-interpolation
	void computeSolutionQI(float nsmooth);

	//	OpenMP version
	// evaluate
	float valueMP(float x, float y, float z, int *tb);
};
