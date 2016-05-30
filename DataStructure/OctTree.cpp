// OctTree.cpp: OctTree 
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
#include <afxwin.h>         //  
#include <afxext.h>         // 
#include <afxdisp.h>        // 
#include <afxdtctl.h>		// 
#ifndef _AFX_NO_AFXCMN_SUPPORT
#include <afxcmn.h>			//
#endif // _AFX_NO_AFXCMN_SUPPORT

#include "windows.h"
#include "OctTree.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#include "math.h"

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////

OctTree::OctTree()
{
	root = NULL;
	table = NULL;
	point_list = NULL;
	normal_list = NULL;
}

OctTree::~OctTree()
{
	if(root != NULL)
		delete root;
	if(table != NULL)
		delete table;
	if(point_list != NULL)
		delete[] point_list;
	if(normal_list != NULL)
		delete[] normal_list;
}

OctTree::Cell::Cell(int level, float ox, float oy, float oz){
	isLeaf = true;
	this->ox = ox;
	this->oy = oy;
	this->oz = oz;
	index_N = 0;
	this->level = level;
}

void OctTree::Cell::createChild(float sizeX, float sizeY, float sizeZ){
	isLeaf = false;
	float sx = 0.5f*getSizeX(sizeX);
	float sy = 0.5f*getSizeY(sizeY);
	float sz = 0.5f*getSizeZ(sizeZ);
	float ox1, oy1, oz1;
	for(int i=0; i<8; i++){
		if(i%2>0)
			ox1 = ox + sx;
		else
			ox1 = ox;

		if(i%4>1)
			oy1 = oy + sy;
		else
			oy1 = oy;

		if(i>3)
			oz1 = oz + sz;
		else
			oz1 = oz;

		child[i] = new Cell(level+1, ox1, oy1, oz1);
	}
}

OctTree::Cell::~Cell(){
	if(!isLeaf){
		for(int i=0; i<8; i++){
			if(child[i] != NULL)
				delete child[i];
		}
	}
	else{
		if(index_N != 0)
			delete[] index;
	}
}

inline float OctTree::Cell::getSizeX(float sizeX)
{
	return (float)(sizeX*pow(2.0, -level));
}

inline float OctTree::Cell::getSizeY(float sizeY)
{
	return (float)(sizeY*pow(2.0, -level));
}

inline float OctTree::Cell::getSizeZ(float sizeZ)
{
	return (float)(sizeZ*pow(2.0, -level));
}

void OctTree::setPointSet(PointSet *ps)
{
	this->ps = ps;
	int point_N = ps->point_N;
	
	float xmin, xmax, ymin, ymax, zmin, zmax;
	ps->getBound(xmin, xmax, ymin, ymax, zmin, zmax);
	float size = max(xmax-xmin, max(ymax-ymin, zmax-zmin));

	/*
	sizeX = size; //(xmax-xmin);
	sizeY = size; //(ymax-ymin);
	sizeZ = size; //(zmax-zmin);
	*/

	sizeX = xmax-xmin;
	sizeY = ymax-ymin;
	sizeZ = zmax-zmin;

	originX = 0.5f*(xmin+xmax) - 0.5*sizeX; //xmin - 0.1f*(xmax-xmin);
	originY = 0.5f*(ymin+ymax) - 0.5*sizeY; //ymin - 0.1f*(ymax-ymin);
	originZ = 0.5f*(zmin+zmax) - 0.5*sizeZ; //zmin - 0.1f*(zmax-zmin);

	root = new Cell(0, originX, originY, originZ);

	float (*point)[3] = ps->point;
	for(int i=0; i<point_N; i++)
		root->addPoint(i, ps, sizeX, sizeY, sizeZ);
	root->cleanChild();
	table = new int[point_N];
}

void OctTree::Cell::addPoint(int i, PointSet* ps, float sizeX, float sizeY, float sizeZ)
{
	if(!isLeaf){
		int j = 0;
		if(ps->point[i][0] >= 0.5*getSizeX(sizeX) + ox)
			j += 1;
		if(ps->point[i][1] >= 0.5*getSizeY(sizeY) + oy)
			j += 2;
		if(ps->point[i][2] >= 0.5*getSizeZ(sizeZ) + oz)
			j += 4;
		child[j]->addPoint(i, ps, sizeX, sizeY, sizeZ);
		return;
	}

	if(index_N == 0)
		index = new int[MAX];
	if(index_N + 1 > MAX){
		createChild(sizeX, sizeY, sizeZ);
		index_N = 0;
		for(int j=0; j<MAX; j++){
			addPoint(index[j], ps, sizeX, sizeY, sizeZ);
		}
		addPoint(i, ps, sizeX, sizeY, sizeZ);
		delete[] index;
		return;
	}
	/*
	for(int j=0; j<index_N; j++){
		float *q = ps->point[index[j]];
		if(p[0] == q[0] && p[1] == q[1] && p[2] == q[2])
			return;
	}
	*/
	index[index_N++] = i;
}

int* OctTree::getIndexTable(int &tableSize, float x, float y, float z, float r)
{
	tableSize = 0;
	root->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
	return table;
}
bool OctTree::getIndexTable(int *table, int &tableSize, float x, float y, float z, float r)
{
	tableSize = 0;
	root->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
	return true;
}
void OctTree::Cell::addIntoTable(PointSet* ps, int *table, int &tableSize, float x, float y, float z, float r, 
								 float sizeX, float sizeY, float sizeZ)
{
	if(isLeaf){
		for(int i=0; i<index_N; i++){
			float *p = ps->point[index[i]];
			if(r*r > (p[0]-x)*(p[0]-x)+(p[1]-y)*(p[1]-y)+(p[2]-z)*(p[2]-z)){
				table[tableSize++] = index[i];
			}
		}
		return;
	}
	float cx = ox + 0.5f*getSizeX(sizeX);
	float cy = oy + 0.5f*getSizeY(sizeY);
	float cz = oz + 0.5f*getSizeZ(sizeZ);
	if(x + r < cx){
		if(y + r < cy){
			if(z + r < cz){
				if(child[0] != NULL)
					child[0]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else if(z - r >= cz){
				if(child[4] != NULL)
					child[4]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else{
				if(child[0] != NULL)
					child[0]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[4] != NULL)
					child[4]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
		}
		else if(y - r >= cy){
			if(z + r < cz){
				if(child[2] != NULL)
					child[2]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else if(z - r >= cz){
				if(child[6] != NULL)
					child[6]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else{
				if(child[2] != NULL)
					child[2]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[6] != NULL)
					child[6]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
		}
		else{
			if(z + r < cz){
				if(child[0] != NULL)
					child[0]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[2] != NULL)
					child[2]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else if(z - r >= cz){
				if(child[4] != NULL)
					child[4]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[6] != NULL)
					child[6]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else{
				if(child[0] != NULL)
					child[0]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[4] != NULL)
					child[4]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[2] != NULL)
					child[2]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[6] != NULL)
					child[6]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
		}
	}
	else if(x -r >= cx){
		if(y + r < cy){
			if(z + r < cz){
				if(child[1] != NULL)
					child[1]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else if(z - r >= cz){
				if(child[5] != NULL)
					child[5]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else{
				if(child[1] != NULL)
					child[1]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[5] != NULL)
					child[5]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
		}
		else if(y - r >= cy){
			if(z + r < cz){
				if(child[3] != NULL)
					child[3]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else if(z - r >= cz){
				if(child[7] != NULL)
					child[7]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else{
				if(child[3] != NULL)
					child[3]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[7] != NULL)
					child[7]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
		}
		else{
			if(z + r < cz){
				if(child[1] != NULL)
					child[1]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[3] != NULL)
					child[3]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else if(z - r >= cz){
				if(child[5] != NULL)
					child[5]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[7] != NULL)
					child[7]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else{
				if(child[1] != NULL)
					child[1]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[5] != NULL)
					child[5]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[3] != NULL)
					child[3]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[7] != NULL)
					child[7]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
		}
	}
	else{
		if(y + r < cy){
			if(z + r < cz){
				if(child[0] != NULL)
					child[0]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[1] != NULL)
					child[1]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else if(z - r >= cz){
				if(child[4] != NULL)
					child[4]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[5] != NULL)
					child[5]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else{
				if(child[0] != NULL)
					child[0]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[1] != NULL)
					child[1]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[4] != NULL)
					child[4]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[5] != NULL)
					child[5]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
		}
		else if(y - r >= cy){
			if(z + r < cz){
				if(child[2] != NULL)
					child[2]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[3] != NULL)
					child[3]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else if(z - r >= cz){
				if(child[6] != NULL)
					child[6]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[7] != NULL)
					child[7]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else{
				if(child[2] != NULL)
					child[2]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[3] != NULL)
					child[3]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[6] != NULL)
					child[6]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[7] != NULL)
					child[7]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
		}
		else{
			if(z + r < cz){
				if(child[0] != NULL)
					child[0]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[1] != NULL)
					child[1]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[2] != NULL)
					child[2]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[3] != NULL)
					child[3]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else if(z - r >= cz){
				if(child[4] != NULL)
					child[4]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[5] != NULL)
					child[5]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[6] != NULL)
					child[6]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[7] != NULL)
					child[7]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
			else{
				if(child[0] != NULL)
					child[0]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[1] != NULL)
					child[1]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[2] != NULL)
					child[2]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[3] != NULL)
					child[3]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[4] != NULL)
					child[4]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[5] != NULL)
					child[5]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[6] != NULL)
					child[6]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
				if(child[7] != NULL)
					child[7]->addIntoTable(ps, table, tableSize, x, y, z, r, sizeX, sizeY, sizeZ);
			}
		}
	}
}

float OctTree::getAverageLeafSize()
{
	int count = 0;
	float sum = 0;

	root->addLeafSize(sizeX, sizeY, sizeZ, count, sum);

	return sum/count;
}

void OctTree::Cell::addLeafSize(float sizeX, float sizeY, float sizeZ, int &count, float &sum)
{
	if(isLeaf){
		if(index_N == 0)
			return;
		count += 1; //index_N;
		//sum += index_N*(getSizeX(sizeX) + getSizeY(sizeY) + getSizeZ(sizeZ))/3.0f;
		float sx = getSizeX(sizeX);
		float sy = getSizeY(sizeY);
		float sz = getSizeZ(sizeZ);
		sum += sqrt(sx*sx+sy*sy+sz*sz);
	}
	else{
		for(int i=0; i<8; i++){
			if(child[i] != NULL)
				child[i]->addLeafSize(sizeX, sizeY, sizeZ, count, sum);
		}
	}
}

void OctTree::computeCenter()
{
	root->computeCenter(ps);
	root->normalizeNormals();
}

void OctTree::Cell::computeCenter(PointSet* ps)
{
	p[0] = p[1] = p[2] = 0;
	n[0] = n[1] = n[2] = 0;
	if(isLeaf){	
		for(int i=0; i<index_N; i++){
			float *q = ps->point[index[i]];
			float *m = ps->normal[index[i]];

			p[0] += q[0];
			p[1] += q[1];
			p[2] += q[2];

			n[0] += m[0];
			n[1] += m[1];
			n[2] += m[2];
		}
	}
	else{

		int i;
		for(i=0; i<8; i++){
			if(child[i] != NULL)
				child[i]->computeCenter(ps);
		}

		index_N = 0;
		for(i=0; i<8; i++){
			if(child[i] == NULL)
				continue;
			int w = child[i]->index_N;
			float *q = child[i]->p;
			float *m = child[i]->n;

			p[0] += w*q[0];
			p[1] += w*q[1];
			p[2] += w*q[2];

			n[0] += w*m[0];
			n[1] += w*m[1];
			n[2] += w*m[2];

			index_N += w;
		}
	}
	if(index_N == 0)
		return;
	p[0] /= index_N;
	p[1] /= index_N;
	p[2] /= index_N;

	n[0] /= index_N;
	n[1] /= index_N;
	n[2] /= index_N;
}

void OctTree::Cell::normalizeNormals()
{
	double len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
	if((float)len != 0){
		n[0] = (float)(n[0]/len);
		n[1] = (float)(n[1]/len);
		n[2] = (float)(n[2]/len);
	}
	if(!isLeaf){
		for(int i=0; i<8; i++){
			if(child[i] != NULL)
				child[i]->normalizeNormals();
		}
	}
}



void OctTree::getPointList(int &size, float **&p_list, float **&n_list, 
						  float x, float y, float z, float r, int max)
{
	size = 0;
	if(point_list == NULL)
		point_list = new float*[ps->point_N];
	if(normal_list == NULL)
		normal_list = new float*[ps->point_N];

	root->addPointIntoList(ps, size, point_list, normal_list, x, y, z, r, max);

	p_list = this->point_list;
	n_list = this->normal_list;
}


void OctTree::Cell::addPointIntoList(PointSet* ps, int &size, float **p_list, float **n_list, 
									 float x, float y, float z, float r, int max)
{
	if(level < max && !isLeaf){
		for(int i=0; i<8; i++){
			if(child[i] != NULL)
				child[i]->addPointIntoList(ps, size, p_list, n_list, x, y, z, r, max);
		}
	}
	else if(level == max){
		if(index_N == 0)
			return;
		if(r*r > (p[0]-x)*(p[0]-x)+(p[1]-y)*(p[1]-y)+(p[2]-z)*(p[2]-z)){
			p_list[size] = p;
			n_list[size] = n;
			size++;
		}
	}
	else if(index_N != 0){
		for(int i=0; i<index_N; i++){
			p_list[size] = ps->point[index[i]];
			n_list[size] = ps->normal[index[i]];
			size++;
		}		
	}
}

void OctTree::Cell::computeCellPoint(PointSet *ps, double *&Q)
{
	Q = new double[10];
	MAT_INIT(Q);
	if(isLeaf){
		if(index_N == 0)
			return;

		int i;
		for(i=0; i<index_N; i++){
			float *q = ps->point[index[i]];
			float *m = ps->normal[index[i]];
			double Qi[10];
			double dot = q[0]*m[0] + q[1]*m[1] + q[2]*m[2];
			double md[3];
			md[0] = m[0];
			md[1] = m[1];
			md[2] = m[2];
			MATRIX(Qi, md, dot);
			MAT_SUM(Q, Qi);
		}
		
		double min = Q_ERR(Q, ps->point[index[0]]);
		int min_i = 0;
		for(i=1; i<index_N; i++){
			double t = Q_ERR(Q, ps->point[index[i]]);
			if(t < min){
				min_i = i;
				min = t;
			}
		}
		float *q = ps->point[index[min_i]];
		p[0] = q[0];
		p[1] = q[1];
		p[2] = q[2];
		float *m = ps->normal[index[min_i]];
		n[0] = m[0];
		n[1] = m[1];
		n[2] = m[2];
	}
	else{
		int i;
		for(i=0; i<8; i++){
			double *Qi;
			child[i]->computeCellPoint(ps, Qi);
			MAT_SUM(Q, Qi);
			delete[] Qi;
		}

		double min = 1000000000;
		int min_i = 0;
		for(i=0; i<8; i++){
			if(child[i] == NULL)
				continue;
			if(child[i]->index_N == 0)
				continue;
			index_N += child[i]->index_N;
			float *q = child[i]->p;
			double t = Q_ERR(Q, q);
			if(t < min){
				min_i = i;
				min = t;
			}
		}
		float *q = child[min_i]->p;
		p[0] = q[0];
		p[1] = q[1];
		p[2] = q[2];
		float *m = child[min_i]->n;
		n[0] = m[0];
		n[1] = m[1];
		n[2] = m[2];
	}
}

void OctTree::computeCellPoint()
{
	double *Q;
	root->computeCellPoint(ps, Q);
	delete[] Q;
}

void OctTree::Cell::splitChild(int M, PointSet *ps, float sizeX, float sizeY, float sizeZ)
{
	if(!isLeaf){
		for(int i=0; i<8; i++){
			if(child[i] != NULL)
				child[i]->splitChild(M, ps, sizeX, sizeY, sizeZ);
		}
	}
	else if(index_N > 1 && level < M-1){
		createChild(sizeX, sizeY, sizeZ);
		int i;
		for(i=0; i<index_N; i++)
			addPoint(index[i], ps, sizeX, sizeY, sizeZ);
		delete[] index;
		for(i=0; i<8; i++){
			if(child[i]->index_N == 0){
				delete child[i];
				child[i] = NULL;
			}
			else
				child[i]->splitChild(M, ps, sizeX, sizeY, sizeZ);
		}
	}
}

void OctTree::splitChild(int M)
{
	root->splitChild(M, ps, sizeX, sizeY, sizeZ);
}

void OctTree::Cell::cleanChild()
{
	if(!isLeaf){
		for(int i=0; i<8; i++){
			if(child[i]->isLeaf && child[i]->index_N == 0){
				delete child[i];
				child[i] = NULL;
			}
			else
				child[i]->cleanChild();
		}
	}
}

/*
void OctTree::setDoubledPoints(int &size, float (*&p_list)[3], float *&value, float off)
{
	p_list = new float[2*ps->point_N][3];
	value = new float[2*ps->point_N];
}

void OctTree::Cell::addDoubledPointIntoList(PointSet *ps, int &size, float (*&p_list )[3], float *&value, float off)
{
	if(!isLeaf){
		for(int i=0; i<8; i++){
			if(child[i] != NULL)
				child[i]->addDoubledPointIntoList(ps, size, p_list, value, off);
		}
	}
	else if(index_N != 0){
		for(int i=0; i<index_N; i++){
			
			p_list[size] = ps->point[index[i]];
			n_list[size] = ps->normal[index[i]];
			size++;
		}
		
	}
}
*/