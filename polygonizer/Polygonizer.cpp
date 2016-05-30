// Polygonizer.cpp: Polygonizer
//
//////////////////////////////////////////////////////////////////////

#include "Polygonizer.h"
#include "../numericalC/SVD.h"
#include "../DataStructure/PolygonalMesh.h"


#include "time.h"
#include "stdio.h"
#include "stdlib.h"
#include "memory.h"
#include "omp.h"
//////////////////////////////////////////////////////////////////////

Polygonizer::Polygonizer()
{
	thread_num = 8;
}

Polygonizer::~Polygonizer()
{

}

void Polygonizer::setDim(int dimX, int dimY, int dimZ)
{
	this->dimX = dimX;
	this->dimY = dimY;
	this->dimZ = dimZ;
}

void Polygonizer::setOrigin(float x, float y, float z)
{
	originX = x;
	originY = y;
	originZ = z;
}

void Polygonizer::setSpace(float x, float y, float z)
{
	spaceX = x;
	spaceY = y;
	spaceZ = z;
}

inline void Polygonizer::weightD(float g[], double d2, double R2, float vx, float vy, float vz)
{
	if(R2 <= d2 || d2 == 0.0){
		g[0] = g[1] = g[2] = 0;
		return;
	}
	double invT2 = 1.0/R2;
	double r = sqrt(d2*invT2);

	double s = 1.0-r;
	double s3 = s*s*s;
	double t = -20.0*s3*invT2;
	g[0] = (float)(vx * t);
	g[1] = (float)(vy * t);
	g[2] = (float)(vz * t);	
}
inline float Polygonizer::value(float vx, float vy, float vz, float nx, float ny, float nz, double R2, bool &isValid)
{
	double d2 = vx*vx+vy*vy+vz*vz;
	if(d2>=R2)
	{
		isValid = false;
		return 0;
	}
	
	isValid = true;
	float g[3];
	weightD(g, d2, R2, vx, vy, vz);
	float f = -(g[0]*nx + g[1]*ny + g[2]*nz);
}

PolygonalMesh* Polygonizer::dualContouring(float epsilon, float tau)
{

	clock_t start = clock();
	//	first pass: remove the grid points outside of the support field

	int (*validID)[3];

//	float p[3];
	bool ***isIn = new bool**[dimZ+1];
	bool ***isValid = new bool**[dimZ+1];
	int ***index = new int**[dimZ+1];

//	float ***value = new float**[dimZ+1];

	int i;
#pragma omp parallel for num_threads(thread_num),schedule(dynamic)
	for(i=0; i<dimZ+1; i++){
		isIn[i] = new bool*[dimY+1];
		isValid[i] = new bool*[dimY+1];
		index[i] = new int*[dimY+1];

//		value[i] = new float*[dimY+1];
		for(int j=0; j<dimY+1; j++){
			isIn[i][j] = new bool[dimX+1];
			isValid[i][j] = new bool[dimX+1];
			index[i][j] = new int[dimX+1];

//			value[i][j] = new float[dimX+1];
			for(int k=0; k<dimX+1; k++)
			{
				isIn[i][j][k] = false;
				isValid[i][j][k] = false;
				index[i][j][k] = -1;

	//			value[i][j][k] = -10000.0f;
			}
		}
	}
	//int num = (dimZ+1)*(dimY+1)*(dimX+1);
	//memset(isIn, 0, num*sizeof(bool));
	//memset(isValid, 0, num*sizeof(bool));

	int count = 0;
#pragma omp parallel for num_threads(thread_num),schedule(dynamic)
	for(i=0; i<dimZ+1; i++)
	{
		float p[3];
		p[2] = originZ + i*spaceZ;
		for(int j=0; j<dimY+1; j++)
		{
			p[1] = originY + j*spaceY;
			for(int k=0; k<dimX+1; k++)
			{
				p[0] = originX + k*spaceX;
				isIn[i][j][k] 
					= (func->value(p[0], p[1], p[2], isValid[i][j][k])>0);

			}
		}
	}

	//clock_t during = clock()-start;
	//printf("Time: %f\n", during*0.001f);

	//printf("The dimension: %d %d %d\n", dimX, dimY, dimZ);
/*	
//---------------------------------------------------
	//	write the distance field to a file

	FILE* file0 = fopen("df_fig.txt","w");
	FILE* file1 = fopen("df1.txt","wb");
	FILE* file2 = fopen("df2.txt","wb");
	FILE* file3 = fopen("df3.txt","wb");
	fprintf(file0, "%f %f %f\n", originX, originY, originZ);
	fprintf(file0, "%f %f %f\n", spaceX, spaceY, spaceZ);
	fprintf(file0, "data format: x_index y_index z_index distanceValue");
	fprintf(file0, "data type: int int int float");
	FILE* file;
	file = file2;
	//fwrite(&originX, sizeof(float), 1, file);	// fprintf(file, "%f %f %f\n", originX, originY, originZ);
	//fwrite(&originY, sizeof(float), 1, file);
	//fwrite(&originZ, sizeof(float), 1, file);

	//fwrite(&spaceX, sizeof(float), 1, file);	//fprintf(file, "%f %f %f", spaceX, spaceY, spaceZ);
	//fwrite(&spaceY, sizeof(float), 1, file);	
	//fwrite(&spaceZ, sizeof(float), 1, file);	
	int tempCount = 0;
	for(i = 0; i< dimZ+1; i++)
	{
		for(int j = 0; j< dimY+1; j++)
		{
			for(int k= 0; k< dimX+1; k++)
			{
				if(value[i][j][k]>-10000.0f)
				{
//					fprintf(file, "\n%d %d %d %f", k, j, i, value[i][j][k]);
					//if(tempCount>=count/3 && tempCount<2*count/3)
					//	file = file2;
					//if(tempCount>=2*count/3)
					//	file = file3;
					fwrite(&k, sizeof(int), 1, file);
					fwrite(&j, sizeof(int), 1, file);
					fwrite(&i, sizeof(int), 1, file);
					fwrite(&(value[i][j][k]), sizeof(float), 1, file);
//					tempCount++;
				}
			}
		}
	}

	fclose(file0);
	fclose(file1);
	fclose(file2);
	fclose(file3);
	////	read the binary file for testing
	//file = fopen("df.txt", "rb");
	//float tempOriginX, tempOriginY, tempOriginZ;
	//float tempSpaceX, tempSpaceY, tempSpaceZ;

	//fread(&tempOriginX, sizeof(float), 1, file);
	//fread(&tempOriginY, sizeof(float), 1, file);
	//fread(&tempOriginZ, sizeof(float), 1, file);
	//fread(&tempSpaceX, sizeof(float), 1, file);
	//fread(&tempSpaceY, sizeof(float), 1, file);
	//fread(&tempSpaceZ, sizeof(float), 1, file);
	//printf("%f %f %f \n %f %f %f\n", tempOriginX, tempOriginY, tempOriginZ, tempSpaceX, tempSpaceY, tempSpaceZ);

	//while(!feof(file))
	//{
	//	if(count > 10)
	//	{
	//		break;
	//	}
	//	int ii[3];
	//	fread(ii, sizeof(int), 3, file);
	//	float tempValue;
	//	fread(&tempValue, sizeof(float), 1, file);
	//	
	//	printf("%d %d %d %f\n", ii[0], ii[1], ii[2], tempValue);
	//	count++;
	//}

	//fclose(file);

	for(i=0; i<dimZ+1; i++){
		for(int j=0; j<dimY+1; j++){
			delete[] value[i][j];
		}
		delete[] value[i];
	}
	delete[] value;
//---------------------------------------------------
*/
	start = clock();
	int current = 0;
	int face_N = 0;
	int face_X_N, face_Y_N, face_Z_N;
	face_X_N = face_Y_N = face_Z_N = 0;
	for(i=0; i<dimZ; i++)
		for(int j=0; j<dimY; j++)
			for(int k=0; k<dimX-1; k++){
				if(!isValid[i][j][k] || !isValid[i][j][k+1])
					continue;
				if((!isIn[i][j][k] && isIn[i][j][k+1]) || 
				    (isIn[i][j][k] && !isIn[i][j][k+1])){
					face_X_N++;
					if(index[i][j][k+1] < 0){
						index[i][j][k+1] = current;
						current++;
					}
					if(index[i][j+1][k+1] < 0){
						index[i][j+1][k+1] = current;
						current++;
					}
					if(index[i+1][j][k+1] < 0){
						index[i+1][j][k+1] = current;
						current++;
					}
					if(index[i+1][j+1][k+1] < 0){
						index[i+1][j+1][k+1] = current;
						current++;
					}
				}
			}

	for(i=0; i<dimX; i++)
		for(int j=0; j<dimZ; j++)
			for(int k=0; k<dimY-1; k++){
				if(!isValid[j][k][i] || !isValid[j][k+1][i])
					continue;
				if((!isIn[j][k][i] && isIn[j][k+1][i]) || 
				    (isIn[j][k][i] && !isIn[j][k+1][i])){
					face_Y_N++;
					if(index[j][k+1][i] < 0){
						index[j][k+1][i] = current;
						current++;
					}
					if(index[j+1][k+1][i] < 0){
						index[j+1][k+1][i] = current;
						current++;
					}
					if(index[j][k+1][i+1] < 0){
						index[j][k+1][i+1] = current;
						current++;
					}
					if(index[j+1][k+1][i+1] < 0){
						index[j+1][k+1][i+1] = current;
						current++;
					}
				}
			}

	for(i=0; i<dimY; i++)
		for(int j=0; j<dimX; j++)
			for(int k=0; k<dimZ-1; k++){
				if(!isValid[k][i][j] || !isValid[k+1][i][j])
					continue;
				if((!isIn[k][i][j] && isIn[k+1][i][j]) || 
				    (isIn[k][i][j] && !isIn[k+1][i][j])){
					face_Z_N++;
					if(index[k+1][i][j] < 0){
						index[k+1][i][j] = current;
						current++;
					}
					if(index[k+1][i+1][j] < 0){
						index[k+1][i+1][j] = current;
						current++;
					}
					if(index[k+1][i][j+1] < 0){
						index[k+1][i][j+1] = current;
						current++;
					}
					if(index[k+1][i+1][j+1] < 0){
						index[k+1][i+1][j+1] = current;
						current++;
					}
				}
			}
	face_N = face_X_N+face_Y_N+face_Z_N;

	//during = 0;
	//during = clock()-start;
	//printf("Time2: %f\n", during*0.001);

	PolygonalMesh* mesh = new PolygonalMesh;
	int vertex_N = current;
	mesh->setVertexCount(vertex_N);
	float (*vertex)[3] = mesh->vertex;
	int* degree = mesh->degree_f = new int[vertex_N];
	current = 0;
	for(i=0; i<vertex_N; i++){
		vertex[i][0] = vertex[i][1] = vertex[i][2] = 0;
		degree[i] = 0;
	}

	mesh->setFaceCount(face_N);
	double (*Q)[10] = new double[vertex_N][10];
	for(i=0; i<vertex_N; i++)
		MAT_INIT(Q[i]);
	for(i=0; i<face_N; i++)
		mesh->setPolygonCount(i, 4);
	int **face = mesh->face;
//	bool flag = false;

	validID = new int[face_N][3];
	face_N = 0;

	//start = clock();
	for(i=0; i<dimZ; i++)
		for(int j=0; j<dimY; j++)
			for(int k=0; k<dimX-1; k++){
				if(!isValid[i][j][k] || !isValid[i][j][k+1])
					continue;
				if(isIn[i][j][k] && !isIn[i][j][k+1]){
					face[face_N][0] = index[i][j][k+1];
					face[face_N][1] = index[i][j+1][k+1];
					face[face_N][2] = index[i+1][j+1][k+1];
					face[face_N][3] = index[i+1][j][k+1];

					validID[face_N][2] = i;	validID[face_N][1] = j;	validID[face_N][0] = k;

					face_N++;
//					flag = true;
				}
				else if(!isIn[i][j][k] && isIn[i][j][k+1]){
					face[face_N][0] = index[i][j][k+1];
					face[face_N][1] = index[i+1][j][k+1];
					face[face_N][2] = index[i+1][j+1][k+1];
					face[face_N][3] = index[i][j+1][k+1];

					validID[face_N][2] = i;	validID[face_N][1] = j;	validID[face_N][0] = k;

					face_N++;
//					flag = true;
				}
				//if(!flag)
				//	continue;
				//flag = false;
			}
	//during = clock()-start;
	//printf("Time2: %f\n", during*0.001);

#pragma omp parallel for num_threads(thread_num)
	for(i = 0; i< face_N; i++)
	{
		float p[3], s[3], e[3];
		int ii, jj, kk;
		ii = validID[i][2];	jj = validID[i][1];	kk = validID[i][0];

		s[0] = originX + kk*spaceX;
		e[0] = originX + (kk+1)*spaceX;
		s[1] = e[1] = originY + jj*spaceY;
		s[2] = e[2] = originZ + ii*spaceZ;
		bisection(p, s, e, epsilon);
		
		//float g[3];
		//func->gradient(g, p[0], p[1], p[2]);
		//double len = PolygonalMesh::LENGTH(g);
		////if((float)len == 0)
		//	//continue;
		//double nor[3];
		//double invLen = 1.0/len;
		//nor[0] = g[0]*invLen;
		//nor[1] = g[1]*invLen;
		//nor[2] = g[2]*invLen;

		//double d = -PolygonalMesh::DOT(nor, p);
		//double Q_tmp[10];
		//MATRIX(Q_tmp, nor, d);

		int i0 = index[ii][jj][kk+1];
		//MAT_SUM(Q[i0], Q_tmp);
		vertex[i0][0] += p[0];
		vertex[i0][1] += p[1];
		vertex[i0][2] += p[2];
		degree[i0]++;

		i0 = index[ii][jj+1][kk+1];
		//MAT_SUM(Q[i0], Q_tmp);
		vertex[i0][0] += p[0];
		vertex[i0][1] += p[1];
		vertex[i0][2] += p[2];
		degree[i0]++;

		i0 = index[ii+1][jj+1][kk+1];
		//MAT_SUM(Q[i0], Q_tmp);
		vertex[i0][0] += p[0];
		vertex[i0][1] += p[1];
		vertex[i0][2] += p[2];
		degree[i0]++;

		i0 = index[ii+1][jj][kk+1];
		//MAT_SUM(Q[i0], Q_tmp);
		vertex[i0][0] += p[0];
		vertex[i0][1] += p[1];
		vertex[i0][2] += p[2];
		degree[i0]++;
	}

	//start = clock();
	for(i=0; i<dimX; i++)
		for(int j=0; j<dimZ; j++)
			for(int k=0; k<dimY-1; k++){
				if(!isValid[j][k][i] || !isValid[j][k+1][i])
					continue;
				if(isIn[j][k][i] && !isIn[j][k+1][i]){
					face[face_N][0] = index[j][k+1][i];
					face[face_N][1] = index[j+1][k+1][i];
					face[face_N][2] = index[j+1][k+1][i+1];
					face[face_N][3] = index[j][k+1][i+1];

					validID[face_N][2] = j;	validID[face_N][1] = k;	validID[face_N][0] = i;

					face_N++;
//					flag = true;
				}
				else if(!isIn[j][k][i] && isIn[j][k+1][i]){
					face[face_N][0] = index[j][k+1][i];
					face[face_N][1] = index[j][k+1][i+1];
					face[face_N][2] = index[j+1][k+1][i+1];
					face[face_N][3] = index[j+1][k+1][i];

					validID[face_N][2] = j;	validID[face_N][1] = k;	validID[face_N][0] = i;

					face_N++;
//					flag = true;
				}
				//if(!flag)
				//	continue;
				//flag = false;
			}
	//during = clock()-start;
	//printf("Time2: %f\n", during*0.001);
#pragma omp parallel for num_threads(thread_num)
	for(i = face_X_N; i< face_N; i++)
	{
		float p[3], s[3], e[3];
		int ii, jj, kk;
		jj = validID[i][2];	kk = validID[i][1];	ii = validID[i][0];
		s[1] = originY + kk*spaceY;
		e[1] = originY + (kk+1)*spaceY;
		s[2] = e[2] = originZ + jj*spaceZ;
		s[0] = e[0] = originX + ii*spaceX;
		bisection(p, s, e, epsilon);
		
		//float g[3];
		//func->gradient(g, p[0], p[1], p[2]);
		//double len = PolygonalMesh::LENGTH(g);
		////if((float)len == 0)
		//	//continue;
		//double nor[3];
		//double invLen = 1.0/len;
		//nor[0] = g[0]*invLen;
		//nor[1] = g[1]*invLen;
		//nor[2] = g[2]*invLen;

		//double d = -PolygonalMesh::DOT(nor, p);
		//double Q_tmp[10];
		//MATRIX(Q_tmp, nor, d);
		
		int i0 = index[jj][kk+1][ii];
		//MAT_SUM(Q[i0], Q_tmp);
		vertex[i0][0] += p[0];
		vertex[i0][1] += p[1];
		vertex[i0][2] += p[2];
		degree[i0]++;

		i0 = index[jj+1][kk+1][ii];
		//MAT_SUM(Q[i0], Q_tmp);
		vertex[i0][0] += p[0];
		vertex[i0][1] += p[1];
		vertex[i0][2] += p[2];
		degree[i0]++;

		i0 = index[jj+1][kk+1][ii+1];
		//MAT_SUM(Q[i0], Q_tmp);
		vertex[i0][0] += p[0];
		vertex[i0][1] += p[1];
		vertex[i0][2] += p[2];
		degree[i0]++;

		i0 = index[jj][kk+1][ii+1];
		//MAT_SUM(Q[i0], Q_tmp);
		vertex[i0][0] += p[0];
		vertex[i0][1] += p[1];
		vertex[i0][2] += p[2];
		degree[i0]++;
	}

	//start = clock();
	for(i=0; i<dimY; i++)
		for(int j=0; j<dimX; j++)
			for(int k=0; k<dimZ-1; k++){
				if(!isValid[k][i][j] || !isValid[k+1][i][j])
					continue;
				if(isIn[k][i][j] && !isIn[k+1][i][j]){
					face[face_N][0] = index[k+1][i][j];
					face[face_N][1] = index[k+1][i][j+1];
					face[face_N][2] = index[k+1][i+1][j+1];
					face[face_N][3] = index[k+1][i+1][j];

					validID[face_N][2] = k;	validID[face_N][1] = i;	validID[face_N][0] = j;

					face_N++;
//					flag = true;
				}
				else if(!isIn[k][i][j] && isIn[k+1][i][j]){
					face[face_N][0] = index[k+1][i][j];
					face[face_N][1] = index[k+1][i+1][j];
					face[face_N][2] = index[k+1][i+1][j+1];
					face[face_N][3] = index[k+1][i][j+1];

					validID[face_N][2] = k;	validID[face_N][1] = i;	validID[face_N][0] = j;

					face_N++;
//					flag = true;
				}
				//if(!flag)
				//	continue;
				//flag = false;
			}

	//during = clock()-start;
	//printf("Time2: %f\n", during*0.001);
#pragma omp parallel for num_threads(thread_num)
	for(i = face_X_N+face_Y_N; i< face_N; i++)
	{
		int ii, jj, kk;
		kk = validID[i][2];	ii = validID[i][1];	jj = validID[i][0];
		float p[3], s[3], e[3];
		s[2] = originZ + kk*spaceZ;
		e[2] = originZ + (kk+1)*spaceZ;
		s[0] = e[0] = originX + jj*spaceX;
		s[1] = e[1] = originY + ii*spaceY;
		bisection(p, s, e, epsilon);
		
		//float g[3];
		//func->gradient(g, p[0], p[1], p[2]);
		//double len = PolygonalMesh::LENGTH(g);
		////if((float)len == 0)
		//	//continue;
		//double nor[3];
		//double invLen = 1.0/len;
		//nor[0] = g[0]*invLen;
		//nor[1] = g[1]*invLen;
		//nor[2] = g[2]*invLen;

		//double d = -PolygonalMesh::DOT(nor, p);
		//double Q_tmp[10];
		//MATRIX(Q_tmp, nor, d);

		int i0 = index[kk+1][ii][jj];
		//MAT_SUM(Q[i0], Q_tmp);
		vertex[i0][0] += p[0];
		vertex[i0][1] += p[1];
		vertex[i0][2] += p[2];
		degree[i0]++;

		i0 = index[kk+1][ii][jj+1];
		//MAT_SUM(Q[i0], Q_tmp);
		vertex[i0][0] += p[0];
		vertex[i0][1] += p[1];
		vertex[i0][2] += p[2];
		degree[i0]++;

		i0 = index[kk+1][ii+1][jj+1];
		//MAT_SUM(Q[i0], Q_tmp);
		vertex[i0][0] += p[0];
		vertex[i0][1] += p[1];
		vertex[i0][2] += p[2];
		degree[i0]++;

		i0 = index[kk+1][ii+1][jj];
		//MAT_SUM(Q[i0], Q_tmp);
		vertex[i0][0] += p[0];
		vertex[i0][1] += p[1];
		vertex[i0][2] += p[2];
		degree[i0]++;
	}
	
	//FOR SVD
			
	//float **A = new float*[4];
	//A[1] = new float[4];
	//A[2] = new float[4];
	//A[3] = new float[4];
	//float *w = new float[4];
	//float **v = new float*[4];
	//v[1] = new float[4];
	//v[2] = new float[4];
	//v[3] = new float[4];
	//float *b = new float[4];
	//float *x = new float[4];
#pragma omp parallel for num_threads(thread_num)
	for(i=0; i<vertex_N; i++){
		if(degree[i] == 0)
			continue;
		double invDeg = 1.0/degree[i];
		vertex[i][0] *= invDeg;
		vertex[i][1] *= invDeg;
		vertex[i][2] *= invDeg;
//continue;
//		A[1][1] = (float)Q[i][0];
//		A[2][1] = A[1][2] = (float)Q[i][1];
//		A[3][1] = A[1][3] = (float)Q[i][2];
//		A[2][2] = (float)Q[i][3];
//		A[2][3] = A[3][2] = (float)Q[i][4];
//		A[3][3] = (float)Q[i][5];
//
//		float Av[3];
//		MAT_BY_VEC(Av, Q[i], vertex[i]);
//		b[1] = -(float)Q[i][6] - Av[0];
//		b[2] = -(float)Q[i][7] - Av[1];
//		b[3] = -(float)Q[i][8] - Av[2];
//
//		SVD::svdcmp(A, 3, 3, w, v);
//		
//		float wmax=0.0f;
//		int k;
//		for (k=1;k<=3;k++)
//			if (fabs(w[k]) > wmax) wmax=(float)fabs(w[k]);
//		if(wmax < 0.0000001)
//			continue;
//		float wmin=wmax*(tau);
//		for (k=1;k<=3;k++){
//			if (fabs(w[k]) < wmin) w[k]=0.0;
//		}
//
//		/*
//		for(int k=1;k<=3;k++)
//			if(fabs(w[k]) < tau) w[k] = 0.0;*/
//			
//
//		SVD::svbksb(A, w, v, 3, 3, b, x);
//		if(fabs(x[1]) > spaceX || fabs(x[2]) > spaceY || fabs(x[3]) > spaceZ)
//			continue;
//
//		mesh->vertex[i][0] += x[1];
//		mesh->vertex[i][1] += x[2];
//		mesh->vertex[i][2] += x[3];
	}

	for(i=0; i<dimZ+1; i++){
		for(int j=0; j<dimY+1; j++){
			delete[] isIn[i][j];
			delete[] index[i][j];
			delete[] isValid[i][j];
		}
		delete[] isIn[i];
		delete[] index[i];
		delete[] isValid[i];
	}
	delete[] isIn;
	delete[] index;
	delete[] Q;
	delete[] isValid;

	delete[] validID;

	//during = clock()-start;
	//printf("Time2: %f\n", during*0.001);

	cutQuad(mesh);

	return mesh;
}


float Polygonizer::smoothGridPointValue(int i, int j, int k, float ***fValue, bool ***pValid)
{
	//	assuming the grid point (i, j, k) is a valid boundary grid point
	float tempV = 0;
	int neiNum = 0;
	int l = i-1;
	//	average the field values of 27 points directly
//#define AVERAGE_27
//#define AVERAGE_6
#define VALID_AVERAGE_27
#ifdef AVERAGE_27
	if(l>=0)
	{
		int h = j-1;
		if(h>=0)
		{
			int m = k-1;
			if(m>=0)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
			m = k;
			if(pValid[l][h][m])
			{
				tempV += fValue[l][h][m];
				neiNum++;
			}
			m = k+1;
			if(m<=dimX)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
		}
		h = j;
		{
			int m = k-1;
			if(m>=0)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
			m = k;
			if(pValid[l][h][m])
			{
				tempV += fValue[l][h][m];
				neiNum++;
			}
			m = k+1;
			if(m<=dimX)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
		}
		h = j+1;
		if(h<=dimY)
		{
			int m = k-1;
			if(m>=0)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
			m = k;
			if(pValid[l][h][m])
			{
				tempV += fValue[l][h][m];
				neiNum++;
			}
			m = k+1;
			if(m<=dimX)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
		}
	}
	l = i;
	{
		int h = j-1;
		if(h>=0)
		{
			int m = k-1;
			if(m>=0)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
			m = k;
			if(pValid[l][h][m])
			{
				tempV += fValue[l][h][m];
				neiNum++;
			}
			m = k+1;
			if(m<=dimX)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
		}
		h = j;
		{
			int m = k-1;
			if(m>=0)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
			m = k;
			if(pValid[l][h][m])
			{
				tempV += fValue[l][h][m];
				neiNum++;
			}
			m = k+1;
			if(m<=dimX)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
		}
		h = j+1;
		if(h<=dimY)
		{
			int m = k-1;
			if(m>=0)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
			m = k;
			if(pValid[l][h][m])
			{
				tempV += fValue[l][h][m];
				neiNum++;
			}
			m = k+1;
			if(m<=dimX)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
		}
	}
	l = i+1;
	if(l<=dimZ)
	{
		int h = j-1;
		if(h>=0)
		{
			int m = k-1;
			if(m>=0)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
			m = k;
			if(pValid[l][h][m])
			{
				tempV += fValue[l][h][m];
				neiNum++;
			}
			m = k+1;
			if(m<=dimX)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
		}
		h = j;
		{
			int m = k-1;
			if(m>=0)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
			m = k;
			if(pValid[l][h][m])
			{
				tempV += fValue[l][h][m];
				neiNum++;
			}
			m = k+1;
			if(m<=dimX)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
		}
		h = j+1;
		if(h<=dimY)
		{
			int m = k-1;
			if(m>=0)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
			m = k;
			if(pValid[l][h][m])
			{
				tempV += fValue[l][h][m];
				neiNum++;
			}
			m = k+1;
			if(m<=dimX)
			{
				if(pValid[l][h][m])
				{
					tempV += fValue[l][h][m];
					neiNum++;
				}
			}
		}
	}
	tempV/=neiNum;
	return tempV;
#endif
	//	average the field values of six points directly
#ifdef AVERAGE_6
	if(l>=0 && pValid[l][j][k])
	{
		tempV += fValue[l][j][k];
		neiNum ++;
	}
	l = j-1;
	if(l>=0 && pValid[i][l][k])
	{
		tempV += fValue[i][l][k];
		neiNum ++;
	}
	l = k-1;
	if(l>=0 && pValid[i][j][l])
	{
		tempV += fValue[i][j][l];
		neiNum ++;
	}

	l = i+1;
	if(l<=dimZ && pValid[l][j][k])
	{
		tempV += fValue[l][j][k];
		neiNum ++;
	}

	l = j+1;
	if(l<=dimY && pValid[i][l][k])
	{
		tempV += fValue[i][l][k];
		neiNum ++;
	}

	l = k+1;
	if(l<=dimX && pValid[i][j][l])
	{
		tempV += fValue[i][j][l];
		neiNum ++;
	}

	if(neiNum>0)
		tempV /= neiNum;
	tempV = (tempV+fValue[i][j][k])*0.5;

	return tempV;
#endif
	//	searching the boundary grid points in the 27 points around a point
#ifdef VALID_AVERAGE_27
	if(fValue[i][j][k]>0)
	{
		if(isBoundaryGridPoint(i-1, j-1, k-1, fValue, pValid) && fValue[i-1][j-1][k-1]>0)
		{
			tempV += fValue[i-1][j-1][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i-1, j-1, k, fValue, pValid) && fValue[i-1][j-1][k]>0)
		{
			tempV += fValue[i-1][j-1][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i-1, j-1, k+1, fValue, pValid) && fValue[i-1][j-1][k+1]>0)
		{
			tempV += fValue[i-1][j-1][k+1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i-1, j, k-1, fValue, pValid) && fValue[i-1][j][k-1]>0)
		{
			tempV += fValue[i-1][j][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i-1, j, k, fValue, pValid) && fValue[i-1][j][k]>0)
		{
			tempV += fValue[i-1][j][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i-1, j, k+1, fValue, pValid) && fValue[i-1][j][k+1]>0)
		{
			tempV += fValue[i-1][j][k+1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i-1, j+1, k-1, fValue, pValid) && fValue[i-1][j+1][k-1]>0)
		{
			tempV += fValue[i-1][j+1][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i-1, j+1, k, fValue, pValid) && fValue[i-1][j+1][k]>0)
		{
			tempV += fValue[i-1][j+1][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i-1, j+1, k+1, fValue, pValid) && fValue[i-1][j+1][k+1]>0)
		{
			tempV += fValue[i-1][j+1][k+1];
			neiNum++;
		}

		if(isBoundaryGridPoint(i, j-1, k-1, fValue, pValid) && fValue[i][j-1][k-1]>0)
		{
			tempV += fValue[i][j-1][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i, j-1, k, fValue, pValid) && fValue[i][j-1][k]>0)
		{
			tempV += fValue[i][j-1][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i, j-1, k+1, fValue, pValid) && fValue[i][j-1][k+1]>0)
		{
			tempV += fValue[i][j-1][k+1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i, j, k-1, fValue, pValid) && fValue[i][j][k-1]>0)
		{
			tempV += fValue[i][j][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i, j, k, fValue, pValid) && fValue[i][j][k]>0)
		{
			tempV += fValue[i][j][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i, j, k+1, fValue, pValid) && fValue[i][j][k+1]>0)
		{
			tempV += fValue[i][j][k+1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i, j+1, k-1, fValue, pValid) && fValue[i][j+1][k-1]>0)
		{
			tempV += fValue[i][j+1][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i, j+1, k, fValue, pValid) && fValue[i][j+1][k]>0)
		{
			tempV += fValue[i][j+1][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i, j+1, k+1, fValue, pValid) && fValue[i][j+1][k+1]>0)
		{
			tempV += fValue[i][j+1][k+1];
			neiNum++;
		}

		if(isBoundaryGridPoint(i+1, j-1, k-1, fValue, pValid) && fValue[i+1][j-1][k-1]>0)
		{
			tempV += fValue[i+1][j-1][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j-1, k, fValue, pValid) && fValue[i+1][j-1][k]>0)
		{
			tempV += fValue[i+1][j-1][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j-1, k+1, fValue, pValid) && fValue[i+1][j-1][k+1]>0)
		{
			tempV += fValue[i+1][j-1][k+1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j, k-1, fValue, pValid) && fValue[i+1][j][k-1]>0)
		{
			tempV += fValue[i+1][j][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j, k, fValue, pValid) && fValue[i+1][j][k]>0)
		{
			tempV += fValue[i+1][j][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j, k+1, fValue, pValid) && fValue[i+1][j][k+1]>0)
		{
			tempV += fValue[i+1][j][k+1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j+1, k-1, fValue, pValid) && fValue[i+1][j+1][k-1]>0)
		{
			tempV += fValue[i+1][j+1][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j+1, k, fValue, pValid) && fValue[i+1][j+1][k]>0)
		{
			tempV += fValue[i][j+1][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j+1, k+1, fValue, pValid) && fValue[i+1][j+1][k+1]>0)
		{
			tempV += fValue[i+1][j+1][k+1];
			neiNum++;
		}
	}
	else
	{
		if(isBoundaryGridPoint(i-1, j-1, k-1, fValue, pValid) && fValue[i-1][j-1][k-1]<=0)
		{
			tempV += fValue[i-1][j-1][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i-1, j-1, k, fValue, pValid) && fValue[i-1][j-1][k]<=0)
		{
			tempV += fValue[i-1][j-1][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i-1, j-1, k+1, fValue, pValid) && fValue[i-1][j-1][k+1]<=0)
		{
			tempV += fValue[i-1][j-1][k+1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i-1, j, k-1, fValue, pValid) && fValue[i-1][j][k-1]<=0)
		{
			tempV += fValue[i-1][j][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i-1, j, k, fValue, pValid) && fValue[i-1][j][k]<=0)
		{
			tempV += fValue[i-1][j][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i-1, j, k+1, fValue, pValid) && fValue[i-1][j][k+1]<=0)
		{
			tempV += fValue[i-1][j][k+1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i-1, j+1, k-1, fValue, pValid) && fValue[i-1][j+1][k-1]<=0)
		{
			tempV += fValue[i-1][j+1][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i-1, j+1, k, fValue, pValid) && fValue[i-1][j+1][k]<=0)
		{
			tempV += fValue[i-1][j+1][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i-1, j+1, k+1, fValue, pValid) && fValue[i-1][j+1][k+1]<=0)
		{
			tempV += fValue[i-1][j+1][k+1];
			neiNum++;
		}

		if(isBoundaryGridPoint(i, j-1, k-1, fValue, pValid) && fValue[i][j-1][k-1]<=0)
		{
			tempV += fValue[i][j-1][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i, j-1, k, fValue, pValid) && fValue[i][j-1][k]<=0)
		{
			tempV += fValue[i][j-1][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i, j-1, k+1, fValue, pValid) && fValue[i][j-1][k+1]<=0)
		{
			tempV += fValue[i][j-1][k+1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i, j, k-1, fValue, pValid) && fValue[i][j][k-1]<=0)
		{
			tempV += fValue[i][j][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i, j, k, fValue, pValid) && fValue[i][j][k]<=0)
		{
			tempV += fValue[i][j][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i, j, k+1, fValue, pValid) && fValue[i][j][k+1]<=0)
		{
			tempV += fValue[i][j][k+1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i, j+1, k-1, fValue, pValid) && fValue[i][j+1][k-1]<=0)
		{
			tempV += fValue[i][j+1][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i, j+1, k, fValue, pValid) && fValue[i][j+1][k]<=0)
		{
			tempV += fValue[i][j+1][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i, j+1, k+1, fValue, pValid) && fValue[i][j+1][k+1]<=0)
		{
			tempV += fValue[i][j+1][k+1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j-1, k-1, fValue, pValid) && fValue[i+1][j-1][k-1]<=0)
		{
			tempV += fValue[i+1][j-1][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j-1, k, fValue, pValid) && fValue[i+1][j-1][k]<=0)
		{
			tempV += fValue[i+1][j-1][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j-1, k+1, fValue, pValid) && fValue[i+1][j-1][k+1]<=0)
		{
			tempV += fValue[i+1][j-1][k+1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j, k-1, fValue, pValid) && fValue[i+1][j][k-1]<=0)
		{
			tempV += fValue[i+1][j][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j, k, fValue, pValid) && fValue[i+1][j][k]<=0)
		{
			tempV += fValue[i+1][j][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j, k+1, fValue, pValid) && fValue[i+1][j][k+1]<=0)
		{
			tempV += fValue[i+1][j][k+1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j+1, k-1, fValue, pValid) && fValue[i+1][j+1][k-1]<=0)
		{
			tempV += fValue[i+1][j+1][k-1];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j+1, k, fValue, pValid) && fValue[i+1][j+1][k]<=0)
		{
			tempV += fValue[i][j+1][k];
			neiNum++;
		}
		if(isBoundaryGridPoint(i+1, j+1, k+1, fValue, pValid) && fValue[i+1][j+1][k+1]<=0)
		{
			tempV += fValue[i+1][j+1][k+1];
			neiNum++;
		}
	}

	tempV /= neiNum;
#endif
	return tempV;
}
bool Polygonizer::isBoundaryGridPoint(int i, int j, int k, float ***fValue, bool ***pValid)
{
	if(i<0 || j<0 || k<0 || i>dimZ || j>dimY || k>dimX)
		return false;

	if(!pValid[i][j][k])
		return false;

	if(fValue[i][j][k]>0)
	{
		int l = i-1;
		if(l>=0 && pValid[l][j][k] && fValue[l][j][k]<=0)
			return true;
		l = i+1;
		if(l<=dimZ && pValid[l][j][k] && fValue[l][j][k]<=0)
			return true;
		l = j-1;
		if(l>=0 && pValid[i][l][k] && fValue[i][l][k]<=0)
			return true;
		l = j+1;
		if(l<=dimY && pValid[i][l][k] && fValue[i][l][k]<=0)
			return true;
		l = k-1;
		if(l>=0 && pValid[i][j][l] && fValue[i][j][l]<=0)
			return true;
		l = k+1;
		if(l<=dimX && pValid[i][j][l] && fValue[i][j][l]<=0)
			return true;
	}
	else
	{
		int l = i-1;
		if(l>=0 && pValid[l][j][k] && fValue[l][j][k]>0)
			return true;
		l = i+1;
		if(l<=dimZ && pValid[l][j][k] && fValue[l][j][k]>0)
			return true;
		l = j-1;
		if(l>=0 && pValid[i][l][k] && fValue[i][l][k]>0)
			return true;
		l = j+1;
		if(l<=dimY && pValid[i][l][k] && fValue[i][l][k]>0)
			return true;
		l = k-1;
		if(l>=0 && pValid[i][j][l] && fValue[i][j][l]>0)
			return true;
		l = k+1;
		if(l<=dimX && pValid[i][j][l] && fValue[i][j][l]>0)
			return true;
	}
	return false;
}
void Polygonizer::bisection(float p[], float start[], float end[], float e)
{
	float p1[3], p2[3], p3[3];
	float f1, f2, f3;
	p1[0] = start[0];
	p1[1] = start[1];
	p1[2] = start[2];
	p2[0] = end[0];
	p2[1] = end[1];
	p2[2] = end[2];
	f1 = func->value(p1[0], p1[1], p1[2]);
	f2 = func->value(p2[0], p2[1], p2[2]);

	//return;
	searchZero(p, p1, p2, f1, f2, e);
	return;

	float edgeLen = PolygonalMesh::DIST(p1, p2);
	for(int j=0; j<10; j++){
		if( edgeLen < 0.001){
			p[0] = p3[0];
			p[1] = p3[1];
			p[2] = p3[2];
			break;
		}

		p3[0] = 0.5f*(p1[0] + p2[0]);
		p3[1] = 0.5f*(p1[1] + p2[1]);
		p3[2] = 0.5f*(p1[2] + p2[2]);

		edgeLen *= 0.5;

		f3 = func->value(p3[0], p3[1], p3[2]);
		if(fabs(f3) < 0.001)//001)
		{
			p[0] = p3[0];
			p[1] = p3[1];
			p[2] = p3[2];
			break;
		}
		else if(f1*f3 >= 0){
			p1[0] = p3[0];
			p1[1] = p3[1];
			p1[2] = p3[2];
			f1 = f3;
		}
		else{
			p2[0] = p3[0];
			p2[1] = p3[1];
			p2[2] = p3[2];
			f2 = f3;
		}
	}
	p[0] = p3[0];
	p[1] = p3[1];
	p[2] = p3[2];

}


void Polygonizer::cutQuad(PolygonalMesh* mesh)
{
	int face_N = mesh->face_N;
	int **face = mesh->face;
	//float (*vertex)[3] = mesh->vertex;
	int **new_face = new int*[2*face_N];
	delete[] mesh->poly_N;
	int *poly_N = mesh->poly_N = new int[2*face_N];
//	int *degree = mesh->degree_f;
	for(int i=0; i<face_N; i++){
		poly_N[2*i] = poly_N[2*i+1] = 3;
		int *nf1 = new_face[2*i] = new int[3];
		int *nf2 = new_face[2*i+1] = new int[3];

		int *f = face[i];
		int i0 = f[0];
		int i1 = f[1];
		int i2 = f[2];
		int i3 = f[3];

		//if(degree[i0] + degree[i2] < degree[i1] + degree[i3]){
		//	nf1[0] = i0;
		//	nf1[1] = i1;
		//	nf1[2] = i2;

		//	nf2[0] = i0;
		//	nf2[1] = i2;
		//	nf2[2] = i3;
		//}
		//else{
			nf1[0] = i0;
			nf1[1] = i1;
			nf1[2] = i3;

			nf2[0] = i1;
			nf2[1] = i2;
			nf2[2] = i3;
//		}

		delete[] f;
	}
	delete[] face;
	mesh->face = new_face;
	mesh->face_N = 2*face_N;
	delete[] mesh->normal_f;
	mesh->normal_f = NULL;
	delete[] mesh->normal;
	mesh->normal = NULL;
}

void Polygonizer::searchZero(float p[], float p1[], float p2[], float f1, float f2, float e)
{
	//regula falsa
	float p3[3], f3;
	if(f1 > 0){
		p3[0] = p1[0];
		p3[1] = p1[1];
		p3[2] = p1[2];

		p1[0] = p2[0];
		p1[1] = p2[1];
		p1[2] = p2[2];

		p2[0] = p3[0];
		p2[1] = p3[1];
		p2[2] = p3[2];

		float swap = f1;
		f1 = f2;
		f2 = swap;
	}

	float dt[3];
	dt[0] = (p2[0] - p1[0]);
	dt[1] = (p2[1] - p1[1]);
	dt[2] = (p2[2] - p1[2]);

	int j;
	double edgeLen2 = dt[0]*dt[0]+dt[1]*dt[1]+dt[2]*dt[2];
	for(j=0; j<5; j++){
		if(edgeLen2 < 0.000001)//0.0001)
		{
			p[0] = p3[0];
			p[1] = p3[1];
			p[2] = p3[2];
			return;
		}

		p3[0] = p1[0] + dt[0]*f1/(f1-f2);
		p3[1] = p1[1] + dt[1]*f1/(f1-f2);
		p3[2] = p1[2] + dt[2]*f1/(f1-f2);
		f3 = func->value(p3[0], p3[1], p3[2]);

		if(fabs(f3) < 0.001)//001)//0.000001)
		{
			p[0] = p3[0];
			p[1] = p3[1];
			p[2] = p3[2];
			return;
		}
		if(f3 < 0.0){
			p1[0] = p3[0];
			p1[1] = p3[1];
			p1[2] = p3[2];

			f1 = f3;
		}
		else{
			p2[0] = p3[0];
			p2[1] = p3[1];
			p2[2] = p3[2];

			f2 = f3;
		}
		dt[0] = (p2[0] - p1[0]);
		dt[1] = (p2[1] - p1[1]);
		dt[2] = (p2[2] - p1[2]);
		edgeLen2 = dt[0]*dt[0]+dt[1]*dt[1]+dt[2]*dt[2];
	}

	//Bisection method
	for(j=0; j<5; j++){
		if(edgeLen2 < 0.000001) //0.000001)
		{
			p[0] = p3[0];
			p[1] = p3[1];
			p[2] = p3[2];
			return;
		}

		p3[0] = 0.5f*(p1[0] + p2[0]);
		p3[1] = 0.5f*(p1[1] + p2[1]);
		p3[2] = 0.5f*(p1[2] + p2[2]);
		f3 = func->value(p3[0], p3[1], p3[2]);
		edgeLen2*=0.5;
		if(fabs(f3) < 0.001)//001)//0.000001)
		{
			p[0] = p3[0];
			p[1] = p3[1];
			p[2] = p3[2];
			return;
		}
		else if(f1*f3 >= 0){
			p1[0] = p3[0];
			p1[1] = p3[1];
			p1[2] = p3[2];
			f1 = f3;
		}
		else{
			p2[0] = p3[0];
			p2[1] = p3[1];
			p2[2] = p3[2];
			f2 = f3;
		}
	}

	p[0] = p3[0];
	p[1] = p3[1];
	p[2] = p3[2];
}

