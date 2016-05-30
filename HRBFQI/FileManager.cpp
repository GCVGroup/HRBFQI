// FileManager.cpp: FileManager 
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FileManager.h"

#include "../DataStructure/PolygonalMesh.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
char CAP(char c){
  if (c >= 'a' && c < 'z')
    return c + ('A' - 'a');
  else
    return c;
}
FileManager::FileManager()
{

}

FileManager::~FileManager()
{

}

void FileManager::setFile(FILE *file, CString file_name, CString file_ext)
{
	this->file =  file;
	this->file_name = file_name;
	this->file_ext = file_ext;
}

void FileManager::save(PolygonalMesh *mesh, float oriCt[3], float scale)
{
	if(file_ext == "obj")
	{
		this->saveWaveFront(mesh,oriCt,scale);
	}
}
void FileManager::saveWaveFront(PolygonalMesh *mesh, float oriCt[3], float scale)
{
	float (*vertex)[3] = mesh->vertex;
	int (**face) = mesh->face;
	int vertex_N = mesh->vertex_N;
	int face_N = mesh->face_N;

	int i;
	for(i=0; i<vertex_N; i++){
		fprintf(file, "v ");
		fprintf(file, "%f ", vertex[i][0]/scale+oriCt[0]);
		fprintf(file, "%f ", vertex[i][1]/scale+oriCt[1]);
		fprintf(file, "%f\n", vertex[i][2]/scale+oriCt[2]);
	}

	for(i=0; i<face_N; i++){
		fprintf(file, "f ");
		fprintf(file, "%d ", face[i][0]+1);
		fprintf(file, "%d ", face[i][1]+1);
		fprintf(file, "%d ", face[i][2]+1);
		fprintf(file, "\n");
	}
}
void FileManager::open(PointSet *ps)
{
	//read cloud data file and perform trianlization.

	if(file_ext == "pwn")
	{
		int point_N;
		fscanf(file, "%d", &point_N);
		ps->setPointSize(point_N);
		float x, y, z;

		int i;
		for(i=0; i<point_N; i++){
			fscanf(file, "%f %f %f", &x, &y, &z);
			ps->setPoint(i, x, y, z);
		}
		for(i=0; i<point_N; i++){
			fscanf(file, "%f %f %f", &x, &y, &z);
			ps->setNormal(i, -x, -y, -z);
		}
	}
}

