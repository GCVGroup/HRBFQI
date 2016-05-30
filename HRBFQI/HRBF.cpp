#include "stdafx.h"
#include "HRBF.h"
#include "math.h"
#include "../numericalC/PBCG.h"
#include "../numericalC/LU.h"
#include <float.h>
#include <iostream>
#include "omp.h"

using namespace std;

#define PI 3.14159265
#define OPEN_MP

HRBF::HRBF(void)
{
	tree = NULL;
	sol = NULL;

	thread_num = THREADS_NUM;	//	THREADS_NUM is defined in StdAfx.h

}

HRBF::~HRBF(void)
{
	if(tree != NULL)
		delete tree;
	if(sol != NULL)
		delete sol;
}
inline void HRBF::weightD(float *w, float g[], double d2, float vx, float vy, float vz)
{
	if(T2<d2)
	{
		*w = 0.0f;
		g[0] = g[1] = g[2] = 0.0f;
		return;
	}
	if(d2 == 0.0)
	{
		*w = 1.0f;
		g[0] = g[1] = g[2] = 0.0f;
		return;
	}

	double invT2 = 1.0/T2;
	double r = sqrt(d2*invT2);
	double s = 1.0-r;
	double s3 = s*s*s;
	double s4 = s3*s;

	*w = s4*(4.0*r + 1.0);
	double t = -20.0*s3*invT2;
	g[0] = (float)(vx * t);
	g[1] = (float)(vy * t);
	g[2] = (float)(vz * t);

}

inline void HRBF::weightD(float g[], double d2, float vx, float vy, float vz)
{
	if(T2 < d2 || d2 == 0.0){
		g[0] = g[1] = g[2] = 0;
		return;
	}
	double invT2 = 1.0/T2;
	double r = sqrt(d2*invT2);

	double s = 1.0-r;
	double s3 = s*s*s;
	double t = -20.0*s3*invT2;
	g[0] = (float)(vx * t);
	g[1] = (float)(vy * t);
	g[2] = (float)(vz * t);
}
inline void HRBF::weightH(float h[], double d2, float vx, float vy, float vz)
{
	if(T2< d2)
	{
		h[0] = h[1] = h[2] = h[3] = h[4] = h[5] = h[6] = h[7] = h[8] = 0.0;
		return;
	}
	if(d2 == 0.0)
	{
		h[0] = h[4] = h[8] = -20.0/T2;
		h[1] = h[2] = h[3] = h[5] = h[6] = h[7] = 0.0;		
		return;
	}

	double r = sqrt(d2/T2);
	double s = 1.0-r;
	double s2 = s*s;

	double t1 = 20.0*s2/(T2*T2*r);
	double t2 = -r*s*T2;
	double vx2 = vx*vx;
	double vy2 = vy*vy;
	double vz2 = vz*vz;
	h[0] = t1*(3.0*vx2+t2);
	h[1] = t1*3.0*vx*vy;
	h[2] = t1*3.0*vx*vz;
	h[3] = h[1];//t1*3.0*vx*vy;
	h[4] = t1*(3.0*vy2+t2);
	h[5] = t1*3.0*vy*vz;
	h[6] = h[2];//t1*3.0*vx*vz;
	h[7] = h[5];//t1*3.0*vy*vz;
	h[8] = t1*(3.0*vz2+t2);
}
inline void HRBF::weightDH(float g[], float h[], double d2, float vx, float vy, float vz)
{
	if(T2<d2)
	{
		g[0] = g[1] = g[2] = 0.0f;
		h[0] = h[1] = h[2] = h[3] = h[4] = h[5] = h[6] = h[7] = h[8] = 0.0;
		return;
	}
	if(d2 == 0.0)
	{
		g[0] = g[1] = g[2] = 0.0f;
		h[0] = h[4] = h[8] = -20.0/T2;
		h[1] = h[2] = h[3] = h[5] = h[6] = h[7] = 0.0;		

		return;
	}

	double invT2 = 1.0/T2;
	double r = sqrt(d2*invT2);
	double s = 1.0-r;
	double s2 = s*s;
	double s3 = s2*s;

	double t = -20.0*s3*invT2;
	g[0] = (float)(vx * t);
	g[1] = (float)(vy * t);
	g[2] = (float)(vz * t);


	double t1 = 20.0*s2/(T2*T2*r);
	double t2 = -r*s*T2;
	double vx2 = vx*vx;
	double vy2 = vy*vy;
	double vz2 = vz*vz;
	h[0] = t1*(3.0*vx2+t2);
	h[1] = t1*3.0*vx*vy;
	h[2] = t1*3.0*vx*vz;
	h[3] = h[1];//t1*3.0*vx*vy;
	h[4] = t1*(3.0*vy2+t2);
	h[5] = t1*3.0*vy*vz;
	h[6] = h[2];//t1*3.0*vx*vz;
	h[7] = h[5];//t1*3.0*vy*vz;
	h[8] = t1*(3.0*vz2+t2);
}
inline void HRBF::weightDH(float *w, float g[], float h[], double d2, float vx, float vy, float vz)
{
	if(T2<d2)
	{
		*w = 0.0f;
		g[0] = g[1] = g[2] = 0.0f;
		h[0] = h[1] = h[2] = h[3] = h[4] = h[5] = h[6] = h[7] = h[8] = 0.0;
		return;
	}
	if(d2 == 0.0)
	{
		*w = 1.0f;
		g[0] = g[1] = g[2] = 0.0f;
		h[0] = h[4] = h[8] = -20.0/T2;
		h[1] = h[2] = h[3] = h[5] = h[6] = h[7] = 0.0;		
		return;
	}
	double invT2 = 1.0/T2;
	double r = sqrt(d2*invT2);
	double s = 1.0-r;
	double s2 = s*s;
	double s3 = s2*s;
	double s4 = s3*s;
	*w = s4*(4.0*r + 1.0);

	double t = -20.0*s3*invT2;
	g[0] = (float)(vx * t);
	g[1] = (float)(vy * t);
	g[2] = (float)(vz * t);


	double t1 = 20.0*s2/(T2*T2*r);
	double t2 = -r*s*T2;
	double vx2 = vx*vx;
	double vy2 = vy*vy;
	double vz2 = vz*vz;
	h[0] = t1*(3.0*vx2+t2);
	h[1] = t1*3.0*vx*vy;
	h[2] = t1*3.0*vx*vz;
	h[3] = h[1];//t1*3.0*vx*vy;
	h[4] = t1*(3.0*vy2+t2);
	h[5] = t1*3.0*vy*vz;
	h[6] = h[2];//t1*3.0*vx*vz;
	h[7] = h[5];//t1*3.0*vy*vz;
	h[8] = t1*(3.0*vz2+t2);
}


inline double HRBF::weight(double d2)
{
	if(T2 < d2)
		return 0;
	
	double r = sqrt(d2/T2);
	return pow(1.0 - r, 4)*(4.0*r + 1.0);
}
float HRBF::value(float x, float y, float z)
{
	double f = 0;
	
	int size, i;

	int table[4000];
	tree->getIndexTable(table, size, x, y, z, support);
	float phi;
	float g[3];
	float *c, *n;
	float vx, vy, vz;
	for(i=0; i<size; i++)
	{
		int j = table[i];
		c = ps->point[j];

		vx = x - c[0];
		vy = y - c[1];
		vz = z - c[2];

		int jj = 4*j;
		double disT = vx*vx + vy*vy + vz*vz;
		//	quasi-interpolation
		weightD(g, disT, vx, vy, vz);
		f-= g[0]*sol[jj+2] + g[1]*sol[jj+3] + g[2]*sol[jj+4];

	}
	return (float)f;
}
float HRBF::value(float x, float y, float z, bool &isValid)
{
	double f = 0;

	int size, i;
	//int* table = tree->getIndexTable(size, x, y, z, support);//
	int table[4000];
	tree->getIndexTable(table, size, x, y, z, support);
	if(size == 0){
		isValid = false;
		return 0;
	}

	isValid = true;

	float phi;
	float g[3];
	float *c, *n;
	float vx, vy, vz;
	for(i=0; i<size; i++)
	{
		int j = table[i];
		c = ps->point[j];

		vx = x - c[0];
		vy = y - c[1];
		vz = z - c[2];

		int jj = 4*j;
		double disT = vx*vx + vy*vy + vz*vz;
		//	quasi-interpolation
		weightD(g, disT, vx, vy, vz);
		f-= g[0]*sol[jj+2] + g[1]*sol[jj+3] + g[2]*sol[jj+4];
	}


	return (float)f;
}

float HRBF::valueMP(float x, float y, float z, int *tb)
{
	double f = 0;

	int size, i, j;
	int* table = tb;
	if(table == NULL)
		table = tree->getIndexTable(size, x, y, z, support);
	else
		tree->getIndexTable(table, size, x, y, z, support);

	float vx, vy, vz;
	double disT;
	float phi;
	float g[3];
	float *c, *n;
	for(i=0; i<size; i++)
	{
		j = table[i];
		c = ps->point[j];
		n = ps->normal[j];

		vx = x - c[0];
		vy = y - c[1];
		vz = z - c[2];

		disT = vx*vx + vy*vy + vz*vz;
		weightD(&phi, g, disT, vx, vy, vz);
		int jj = 4*j;

		f+= phi*sol[jj+1] - g[0]*sol[jj+2] - g[1]*sol[jj+3] - g[2]*sol[jj+4];
	}


	return (float)f;
}


void HRBF::gradient(float g[], float x, float y, float z)
{
	g[0] = g[1] = g[2] = 0.0f;
	int size;
//	int* table = tree->getIndexTable(size, x, y, z, support);
	int table[4000];
	tree->getIndexTable(table, size, x, y, z, support);

	float vx, vy, vz;
	float *c;
	double g1[3];
	g1[0] = g1[1] = g1[2] = 0;
	for(int i=0; i<size; i++){
		int j = table[i];
		c = ps->point[j];
		vx = x - c[0];
		vy = y - c[1];
		vz = z - c[2];
		double d2 = vx*vx + vy*vy + vz*vz;
		float hw[9];
		weightH(hw, d2, vx, vy, vz);

		int jj = 4*j;
		g1[0] -= sol[jj+2]*hw[0]+sol[jj+3]*hw[1]+sol[jj+4]*hw[2];
		g1[1] -= sol[jj+2]*hw[3]+sol[jj+3]*hw[4]+sol[jj+4]*hw[5];
		g1[2] -= sol[jj+2]*hw[6]+sol[jj+3]*hw[7]+sol[jj+4]*hw[8];
	}
	g[0] += (float)g1[0];
	g[1] += (float)g1[1];
	g[2] += (float)g1[2];		
}
float HRBF::valueAndGradient(float g[], float x, float y, float z)
{
	g[0] = g[1] = g[2] = 0.0f;

	double f = 0;

	int size, j, i;
	int* table = tree->getIndexTable(size, x, y, z, support);
	
	float vx, vy, vz;
	double d2;
	float *c;
	float gw[3], hw[9], w;
	double g1[3];
	g1[0] = g1[1] = g1[2] = 0;
	double f1 = 0.0;

	for(i=0; i<size; i++)
	{
		j = table[i];
		int jj = 4*j;

		c = ps->point[j];

		vx = x - c[0];
		vy = y - c[1];
		vz = z - c[2];

		d2 = vx*vx + vy*vy + vz*vz;
		weightDH(&w, gw, hw, d2, vx, vy, vz);
		
		f1+= w*sol[jj+1] - gw[0]*sol[jj+2] - gw[1]*sol[jj+3] - gw[2]*sol[jj+4];
		g1[0] += sol[jj+1]*gw[0]-sol[jj+2]*hw[0]-sol[jj+3]*hw[1]-sol[jj+4]*hw[2];
		g1[1] += sol[jj+1]*gw[1]-sol[jj+2]*hw[3]-sol[jj+3]*hw[4]-sol[jj+4]*hw[5];
		g1[2] += sol[jj+1]*gw[2]-sol[jj+2]*hw[6]-sol[jj+3]*hw[7]-sol[jj+4]*hw[8];
	}
	f += f1;
	g[0] += (float)g1[0];
	g[1] += (float)g1[1];
	g[2] += (float)g1[2];		

	return (float)f;
}
void HRBF::setSupport(float T)
{
	this->support = T;
	T2 = T*T;
}
//Octree oprtations
void HRBF::setPointSet(PointSet* ps)
{
	this->ps = ps;
	tree = new OctTree();
	tree->setPointSet(ps);
}
float HRBF::getAveragedLeafSize()
{
	return tree->getAverageLeafSize();
}
int HRBF::getMaximalNeighborsInSupport(float T)
{
	int maxNeiNum = 0;
	int i;
#ifdef OPEN_MP
	int *maxNum = new int[thread_num];
	int **neiTable = new int*[thread_num];
	int estNeiNum = 4000;
	for(i = 0; i< thread_num; i++)
	{
		maxNum[i] = 0;
		neiTable[i] = new int[estNeiNum];
	}

#pragma omp parallel for schedule(dynamic)
	for(i = 0; i< ps->point_N; i++)
	{
		int size = 0;
		float *c = ps->point[i];
		 
		int thread_index = omp_get_thread_num();

		tree->getIndexTable(neiTable[thread_index], size, c[0], c[1], c[2], T);

		if(maxNum[thread_index] < size)
			maxNum[thread_index] = size;
	}
	for(i = 0; i< thread_num; i++)
	{
		if(maxNeiNum<maxNum[i])
			maxNeiNum = maxNum[i];

		delete[] neiTable[i];
	}
	delete[] maxNum;
	delete[] neiTable;
#else
	{
		for(i = 0; i< ps->point_N; i++)
		{
			int size = 0;
			float *c = ps->point[i];
			int* table = tree->getIndexTable(size, c[0], c[1], c[2], T);
			if(maxNeiNum < size)
				maxNeiNum = size;
		}
	}
#endif
	return maxNeiNum;
}

//Fitting

void HRBF::fit(float support, float nsmooth)
{
	this->support = support;
	T2 = support*support;
	this->computeSolutionQI(nsmooth);
}
	
void HRBF::computeSolution(float nsmooth)
{
	return;

	int i;
	int point_N = ps->point_N;

	int coefNum = point_N*4;//+4;
	if(sol != NULL)
		delete[] sol;
	sol = new double[coefNum+1];

	float (*point)[3] = ps->point;
	float (*normal)[3] = ps->normal;
	int total = 0;

	for(i=0; i<point_N; i++){
		int size = 0;
		float *c = point[i];
		int* t = tree->getIndexTable(size, c[0], c[1], c[2], support);
		total += size;
	}
	
	unsigned long totalNum = total*16 - point_N*12;// + point_N*14;
	cout<<"totalNum = "<< totalNum<<endl;
	PBCG solver;
	double *sa = solver.sa = new double[totalNum+2];
	unsigned long *ija = solver.ija = new unsigned long[totalNum+2];
	double *b = new double[coefNum+1];

	ija[1]=coefNum+2;
	unsigned long kk=coefNum+1;

//	int nN = 59;	//for dinasour: 30 61 112; red-circle:13 33 59; armadillo: 33 62 107; buste:33 59 106;dancer:25  62 125;
	double tempNSmooth = nsmooth; //10.0*nN*(1.25/support+483.0/(25.0*T2))*nsmooth;//1000.0/T2*nsmooth;//20000.0*nsmooth;5*(1.0+20.0/(T2));

	////------------------------------
	////	for outputing the matrix "sa"
	//double **sa_bk = new double* [coefNum];
	//for(i = 0; i< coefNum; i++)
	//{
	//	sa_bk[i] = new double [coefNum];
	//	memset(sa_bk[i], 0, sizeof(double)*coefNum);
	//}
	//
	////------------------------------
	for(i=0; i<point_N; i++){
		int ii = 4*i;
		float *n = ps->normal[i];
		double f_sub = 0;
		float g[3] = {0,0,0};
		double gLen2 = g[0]*g[0]+g[1]*g[1]+g[2]*g[2];
		//double gLen = sqrt(gLen2);
		//if(gLen == 0.0)
		//	gLen = 1.0;
		//double dotNorm = n[0]*g[0] + n[1]*g[1] + n[2]*g[2];
		//if(dotNorm<0.0)
		//{
		//	g[0]*=-1.0f;
		//	g[1]*=-1.0f;
		//	g[2]*=-1.0f;
		//}
		if(gLen2 == 0.0)
			g[0] = g[1] = g[2] = 0.0f;
		else
		{
			double invGLen = 1.0/sqrt(gLen2);
			g[0]*=invGLen;
			g[1]*=invGLen;
			g[2]*=invGLen;
		}
		b[ii+1] = -f_sub;
		b[ii+2] = n[0]-g[0];
		b[ii+3] = n[1]-g[1];
		b[ii+4] = n[2]-g[2];
		sol[ii+1] = -f_sub;
		sol[ii+2] = n[0]-g[0];
		sol[ii+3] = n[1]-g[1];
		sol[ii+4] = n[2]-g[2];

		int size = 0;
		float *c = point[i];
		int* table = tree->getIndexTable(size, c[0], c[1], c[2], support);
		int* t = new int[size];
		int m = 0;

		int j;
		for(j=0; j<size; j++){
			t[m++] = table[j];
		}

		for(j=0; j<m; j++){
			for(int k=j+1; k<m; k++){
				if(t[j] > t[k]){
					int tmp = t[j];
					t[j] = t[k];
					t[k] = tmp;
				}
			}
		}
		int pos_x = 4*(m-1);
		int pos_y = 8*(m-1);
		int pos_z = 12*(m-1);
		double sumW = 0;
		float sumWD[3] = {0.0f, 0.0f, 0.0f};
		
		//for(j=0; j<m; j++)
		//{
		//	int k = t[j];
		//	float vx = point[i][0] - point[k][0];
		//	float vy = point[i][1] - point[k][1];
		//	float vz = point[i][2] - point[k][2];
		//	double d2 = vx*vx + vy*vy + vz*vz;
		//	double w = weight(d2);
		//	sumW+=w;
		//	float g[3];
		//	weightD(g, d2, vx, vy, vz);
		//	sumWD[0]+=g[0];
		//	sumWD[1]+=g[1];
		//	sumWD[2]+=g[2];
		//}
		double invSumW = 1.0f;///sumW;
		//double normSumWD[3];
		//normSumWD[0] = sumWD[0]*invSumW;
		//normSumWD[1] = sumWD[1]*invSumW;
		//normSumWD[2] = sumWD[2]*invSumW;

		sa[ii+1] = 1.0 + tempNSmooth;//fsmooth;//(weight(0) + smooth)*invSumW;
		sa[ii+2] = 20.0/T2 + tempNSmooth; //1000.0/T2*nsmooth;//sa[ii+1];
		//
		//
		//
		sa[ii+3] = sa[ii+2];//sa[ii+1];
		sa[ii+4] = sa[ii+2];//sa[ii+1];
		double invSumW2 = 1.0;//invSumW*invSumW;
		////--------------------------------------------------------
		////	write the matrix sa_bk
		//for(j = 0; j<m; j++)
		//{
		//	int h = t[j];
		//	int hh = 4*h;
		//	if(h == i)
		//	{
		//		sa_bk[ii][hh] = 1.0+tempNSmooth;
		//		sa_bk[ii+1][hh+1] = 20.0/T2 + tempNSmooth;
		//		sa_bk[ii+2][hh+2] = sa_bk[ii+1][hh+1];
		//		sa_bk[ii+3][hh+3] = sa_bk[ii+1][hh+1];
		//	}
		//	else
		//	{
		//		float vx = point[i][0] - point[h][0];
		//		float vy = point[i][1] - point[h][1];
		//		float vz = point[i][2] - point[h][2];
		//		double d2 = vx*vx + vy*vy + vz*vz;
		//		double w = weight(d2);
		//		float gw[3];
		//		weightD(gw, d2, vx, vy, vz);
		//		float hw[9];
		//		weightH(hw, d2, vx, vy, vz);
		//		double normW = w;//*invSumW;
//
		//		sa_bk[ii][hh] = w;			sa_bk[ii][hh+1] = -gw[0];		sa_bk[ii][hh+2] = -gw[1];		sa_bk[ii][hh+3] = -gw[2];
		//		sa_bk[ii+1][hh] = gw[0];	sa_bk[ii+1][hh+1] = -hw[0];		sa_bk[ii+1][hh+2] = -hw[1];		sa_bk[ii+1][hh+3] = -hw[2];
		//		sa_bk[ii+2][hh] = gw[1];	sa_bk[ii+2][hh+1] = -hw[3];		sa_bk[ii+2][hh+2] = -hw[4];		sa_bk[ii+2][hh+3] = -hw[5];
		//		sa_bk[ii+3][hh] = gw[2];	sa_bk[ii+3][hh+1] = -hw[6];		sa_bk[ii+3][hh+2] = -hw[7];		sa_bk[ii+3][hh+3] = -hw[8];

		//	}
		//}
		////--------------------------------------------------------

		for(j=0; j<m; j++){
			int k = t[j];
			if(k == i)
			{
				continue;
			}
			float vx = point[i][0] - point[k][0];
			float vy = point[i][1] - point[k][1];
			float vz = point[i][2] - point[k][2];
			double d2 = vx*vx + vy*vy + vz*vz;
			double w = weight(d2);
			float gw[3];
			weightD(gw, d2, vx, vy, vz);
			float hw[9];
			weightH(hw, d2, vx, vy, vz);
			//w *= weightNormal(n, normal[i]);
//			if(w != 0.0)
			double normW = w;//*invSumW;
			int k4 = 4*k;
			{
				//	const part
				kk++;
				sa[kk] = normW;
				ija[kk] = k4+1;

				sa[pos_x+kk] = gw[0];//*invSumW-normW*normSumWD[0];
				ija[pos_x+kk] = k4+1;

				sa[pos_y+kk] = gw[1];//*invSumW-normW*normSumWD[1];
				ija[pos_y+kk] = k4+1;

				sa[pos_z+kk] = gw[2];//*invSumW-normW*normSumWD[2];
				ija[pos_z+kk] = k4+1;

				//	x part
				kk++;
				sa[kk] = -gw[0];//normW*vx;
				ija[kk] = k4+2;

				sa[pos_x+kk] = -hw[0];//normW + vx*g[0]*invSumW - vx*normW*normSumWD[0];
				ija[pos_x+kk] = k4+2;

				sa[pos_y+kk] = -hw[3];//vx*g[1]*invSumW - vx*normW*normSumWD[1];
				ija[pos_y+kk] = k4+2;

				sa[pos_z+kk] = -hw[6];//vx*g[2]*invSumW - vz*normW*normSumWD[2];
				ija[pos_z+kk] = k4+2;

				//	y part
				kk++;
				sa[kk] = -gw[1];//normW*vy;
				ija[kk] = k4+3;

				sa[pos_x+kk] = -hw[1];//vy*g[0]*invSumW - vy*normW*normSumWD[0];
				ija[pos_x+kk] = k4+3;

				sa[pos_y+kk] = -hw[4];//normW + vy*g[1]*invSumW - vy*normW*normSumWD[1];
				ija[pos_y+kk] = k4+3;

				sa[pos_z+kk] = -hw[7];//vy*g[2]*invSumW - vy*normW*normSumWD[2];
				ija[pos_z+kk] = k4+3;
				
				//	z part
				kk++;
				sa[kk] = -gw[2];//normW*vz;
				ija[kk] = k4+4;

				sa[pos_x+kk] = -hw[2];//vz*g[0]*invSumW - vz*normW*normSumWD[0];
				ija[pos_x+kk] = k4+4;

				sa[pos_y+kk] = -hw[5];//vz*g[1]*invSumW - vz*normW*normSumWD[1];
				ija[pos_y+kk] = k4+4;

				sa[pos_z+kk] = -hw[8];//normW + vz*g[2]*invSumW - vz*normW*normSumWD[2];
				ija[pos_z+kk] = k4+4;

			}
		}

		ija[ii+2]=kk+1;
		ija[ii+3]=pos_x+kk+1;
		ija[ii+4]=pos_y+kk+1;
		ija[ii+5]=pos_z+kk+1;
		
		kk = pos_z+kk;

		if(m != 0)
			delete[] t;

	}
	cout<<"kk = "<< kk<<endl;

	int iter;
	double err;
	solver.linbcg(coefNum, b, sol, 4, 1e-6, 1000, &iter, &err);

	delete[] b;
//	//-------------------------------------------------------------------------------------------
//	//	output the matrix sa_bk into a file
//
//	FILE* file = fopen("X:\\Temp\\November\\HermiteRBF\\Results\\20150820\\matrix_ida_r0.8.txt","w");
//	for(i=0; i<coefNum; i++)
//	{
////		fprintf(file, "v ");
//		int j;
//		for(j = 0; j< coefNum; j++)
//		{
//			fprintf(file, "%f ", sa_bk[i][j]);
//		}
//		fprintf(file, "\n");
//	}	
//	fclose(file);
//
//	for(i = 0; i< coefNum; i++)
//		delete[] sa_bk[i];
//	delete[] sa_bk;
//	//-------------------------------------------------------------------------------------------
}




//	quasi-interpolation
void HRBF::computeSolutionQI(float nsmooth)
{
	int i;
	int point_N = ps->point_N;

	int coefNum = point_N*4;
	if(sol != NULL)
		delete[] sol;
	sol = new double[coefNum+1];

	float (*point)[3] = ps->point;
	float (*normal)[3] = ps->normal;

	double regular_coef = 1.0;///(nsmooth + 20.0/T2);
#pragma omp parallel for schedule(dynamic)
	for(i = 0; i< point_N; i++)
	{
		int ii = 4*i;
		sol[ii+1] = 0.0;
		sol[ii+2] = normal[i][0]*regular_coef;
		sol[ii+3] = normal[i][1]*regular_coef;
		sol[ii+4] = normal[i][2]*regular_coef;
	}

}



//2D mesh for cross sections. (for the paper)
PolygonalMesh* HRBF::generateCrossSection(float o[], float t1[], float t2[], int n, int m)
{
	PolygonalMesh *mesh = new PolygonalMesh;
	mesh->setVertexCount(n*m);
	mesh->setFaceCount((n-1)*(m-1));
	int **face = mesh->face;
	float (*vertex)[3] = mesh->vertex;
	float *v = mesh->value = new float[mesh->vertex_N];
	mesh->isValid = new bool[mesh->vertex_N];
	mesh->isCovered = new bool[mesh->vertex_N];

	double t3[3];
	PolygonalMesh::CROSS(t3, t1, t2);
	double len = PolygonalMesh::LENGTH(t3);
	//0.3 for kernal
	t3[0] /= -0.4*len;
	t3[1] /= -0.4*len;
	t3[2] /= -0.4*len;

	int i;
	for(i = 0; i< mesh->vertex_N; i++)
		mesh->isValid[i] = true;

	for(i=0; i<n; i++){
		float p1[3];
		p1[0] = o[0] + t1[0]*i;
		p1[1] = o[1] + t1[1]*i;
		p1[2] = o[2] + t1[2]*i;
		for(int j=0; j<m; j++){
			float p2[3];
			p2[0] = t2[0]*j;
			p2[1] = t2[1]*j;
			p2[2] = t2[2]*j;

			int index = i*m+j;
			vertex[index][0] = p1[0] + p2[0];
			vertex[index][1] = p1[1] + p2[1];
			vertex[index][2] = p1[2] + p2[2];

			bool flag;
			float f = value(vertex[index][0], vertex[index][1], vertex[index][2], flag);//);//
			

			if(_isnan(f)){
				v[index] = 0;
				continue;
			}

			/*
			vertex[index][0] += f*(float)t3[0];
			vertex[index][1] += f*(float)t3[1];
			vertex[index][2] += f*(float)t3[2];
			*/
			
			v[index] = f;
			mesh->isCovered[index] = flag;
			//if(!flag)
				//v[index] = 1000000;
			//else
				//v[index] = -10000000;
		}
	}

	for(i=0; i<n-1; i++){
		for(int j=0; j<m-1; j++){
			int index = i*(m-1)+j;
			mesh->setPolygonCount(index, 4);
			face[index][3] = i*m+j;
			face[index][2] = (i+1)*m+j;
			face[index][1] = (i+1)*m+j+1;
			face[index][0] = i*m+j+1;
			//bool flag = (v[face[index][0]] > 0);
			
			/*
			if(PolygonalMesh::LENGTH(vertex[face[index][0]]) > 5 ||
				PolygonalMesh::LENGTH(vertex[face[index][1]]) > 5 ||
				PolygonalMesh::LENGTH(vertex[face[index][2]]) > 5 ||
				PolygonalMesh::LENGTH(vertex[face[index][3]]) > 5)
				face[index][0] = -1;
				*/
			/*
			for(int k=1; k<4; k++){ 
				if(flag != (v[face[index][k]] > 0)){
					v[i*n+j] = -100;
				}
			}*/
			/*
			if(v[face[index][0]] == v[face[index][1]] 
				== v[face[index][2]] == v[face[index][3]] == 0)
				v[i*n+j] = 1000000;*/
		}
	}
	mesh->computeFaceNormal();
	return mesh;
}

void HRBF::differenceFuncGradientAndNormalPt(double &maxAng, double &aveAng)
{
	double *grad = new double[3*ps->point_N];
	memset(grad, 0, 3*ps->point_N*sizeof(double));

	int i;
	int numTotal = 0;
	double sumE = 0.0;
	double maxE = 0.0;
	for(i = 0; i< ps->point_N; i++)
	{
		int kk = 3*i;
		float g[3];
		float *p = ps->point[i];
		float *n = ps->normal[i];
		gradient(g, p[0], p[1], p[2]);
		double len = PolygonalMesh::LENGTH(g);
		if(len <= 1.0e-10)
			continue;
		double invLen = 1.0/len;
		grad[kk] = g[0]*invLen;
		grad[kk+1] = g[1]*invLen;
		grad[kk+2] = g[2]*invLen;

		double normalDot = (grad[kk]*n[0]+grad[kk+1]*n[1]+grad[kk+2]*n[2]);
		double angle = fabs(acos(normalDot)*180.0/PI);
		if(maxE < angle)
			maxE = angle;
		if(angle > 0.0 && angle <= 180.0)
		{
			sumE += angle;
			numTotal += 1;
		}
	}
	maxAng = maxE;
	aveAng = sumE /(double)numTotal;
	delete[] grad;
}