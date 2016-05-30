// PolygonalMesh.cpp: PolygonalMesh 
//
//////////////////////////////////////////////////////////////////////
#include <afxwin.h>         
#include <afxext.h>         
#include <afxdisp.h>        
#include <afxdtctl.h>		
#ifndef _AFX_NO_AFXCMN_SUPPORT
#include <afxcmn.h>			
#endif // _AFX_NO_AFXCMN_SUPPORT

#include "PolygonalMesh.h"
#include "..\PQP\include\PQP.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////

PolygonalMesh::PolygonalMesh()
{
	face_N = 0;
	vertex_N = 0;

	vertex = NULL;
	face = NULL;

	degree_f = NULL;

	isBound = NULL;

	normal_f = NULL;
	normal = NULL;

	center = NULL;

	value = NULL;
	isValid = NULL;
	isCovered = NULL;
}

PolygonalMesh::~PolygonalMesh()
{
	if(vertex != NULL)
		delete[] vertex;
	if(face != NULL){
		for(int i=0; i<face_N; i++){
			if(poly_N != 0)
				delete[] face[i];
		}
		delete[] face;
	}
	delete[] poly_N;
	if(normal != NULL)
		delete[] normal;
	if(normal_f != NULL)
		delete[] normal_f;

	if(center != NULL)
		delete[] center;

	if(value != NULL)
		delete[] value;
	if(isValid != NULL)
		delete[] isValid;
	if(isCovered != NULL)
		delete[] isCovered;

}

void PolygonalMesh::setupTopologyInfo()	//just for triangular mesh
{
	int i;
	m_vertLinkInfo.resize(vertex_N);
	m_faceLinkInfo.resize(face_N);

	//	get the triangles linked to the vertices
	for(i = 0; i< face_N; i++)
	{
		int j;
		for(j = 0; j< 3; j++)
		{
			m_vertLinkInfo[face[i][j]].m_linkFaceIds.push_back(i);
		}
		m_faceLinkInfo[i].m_isBound = false;
		m_faceLinkInfo[i].m_linkEdgeIds.resize(3);
	}
	
	//	setup the edge link information according the vertices link information
	for(i = 0; i< vertex_N; i++)
	{
		m_vertLinkInfo[i].m_isBound = false;
		INTVECTOR faceIds = m_vertLinkInfo[i].m_linkFaceIds;
		int j;
		for(j = 0; j< faceIds.size(); j++)
		{
			int faceId = faceIds[j];

			int k, vid2;
			for(k = 0; k< 3; k++)
			{
				if(i!=face[faceId][k])	continue;
				vid2 = face[faceId][(k+1)%3];
				if(vid2 > i)	// not record before
				{
					EDGE_LINK_INFO edge_link_info;
					edge_link_info.m_linkVertIds.push_back(i);
					edge_link_info.m_linkVertIds.push_back(vid2);
					edge_link_info.m_linkFaceIds.push_back(faceId);
					int jj;
					for(jj = 0; jj< faceIds.size(); jj++)
					{
						if(jj == j)	continue;	// not same triangle
						int faceId2 = faceIds[jj];
						int kk;
						for(kk = 0; kk< 3; kk++)
						{
							if(vid2 == face[faceId2][kk])
							{
								edge_link_info.m_linkFaceIds.push_back(faceId2);
								break;
							}
						}
						if(edge_link_info.m_linkFaceIds.size()==2)
							break;
					}
					if(edge_link_info.m_linkFaceIds.size()==1)
						edge_link_info.m_isBound = true;
					else
						edge_link_info.m_isBound = false;
					m_edgeLinkInfo.push_back(edge_link_info);
				}
				break;
			}
		}
	}

	//	setup other link information for vertices and faces according the edges information
	for(i = 0; i< m_edgeLinkInfo.size(); i++)
	{
		EDGE_LINK_INFO edge_link_info = m_edgeLinkInfo[i];

		int vid1 = edge_link_info.m_linkVertIds[0];
		int vid2 = edge_link_info.m_linkVertIds[1];

		int j;
		for(j = 0; j< edge_link_info.m_linkVertIds.size(); j++)	// vertices
		{
			int vid = edge_link_info.m_linkVertIds[j];
			m_vertLinkInfo[vid].m_linkEdgeIds.push_back(i);
			
			m_vertLinkInfo[vid].m_isBound = m_vertLinkInfo[vid].m_isBound || edge_link_info.m_isBound;
		}
		
		for(j = 0; j< edge_link_info.m_linkFaceIds.size(); j++)	//	faces
		{
			int fid = edge_link_info.m_linkFaceIds[j];
			int k;
			for(k = 0; k< 3; k++)
			{
				if(face[fid][k] != vid1 && face[fid][k] != vid2)
				{
					m_faceLinkInfo[fid].m_linkEdgeIds[k] = i;
					break;
				}
			}
			m_faceLinkInfo[fid].m_isBound = m_faceLinkInfo[fid].m_isBound || edge_link_info.m_isBound;
		}
	}

	//	setup link face information for face according the face information
//#pragma omp parallel for num_threads(OPENMP_THREADS_NUMBER)
	for(i = 0; i< face_N; i++)
	{
		int j;
		for(j = 0; j< 3; j++)
		{
			int eid = m_faceLinkInfo[i].m_linkEdgeIds[j];
			EDGE_LINK_INFO edge_link_info = m_edgeLinkInfo[eid];

			if(edge_link_info.m_isBound)	continue;

			int fid;
			fid = (i==edge_link_info.m_linkFaceIds[0]?edge_link_info.m_linkFaceIds[1]:edge_link_info.m_linkFaceIds[0]);

			m_faceLinkInfo[i].m_linkFaceIds.push_back(fid);
		}
	}

	//	check the topology
	//	if the model is closed
	//bool isNotEnclosed = false;
	//for(i = 0; i< face_N; i++)
	//{
	//	if(m_faceLinkInfo[i].m_isBound)
	//	{
	//		isNotEnclosed = true;
	//		break;
	//	}
	//}
	//if(isNotEnclosed)
	//{
	//	cout<<"The surface is not enclosed!"<<endl;
	//}
	//else
	//	cout<<"The surface is enclosed!"<<endl;

}

void PolygonalMesh::setVertexCount(int vertex_N)
{
	if(vertex != NULL)
		delete[] vertex;
	this->vertex_N = vertex_N;
	vertex = new float[vertex_N][3];
}

void PolygonalMesh::setFaceCount(int face_N)
{
	if(face != NULL){
		delete[] face;
		delete[] poly_N;
	}
	this->face_N = face_N;
	face = new int*[face_N];
	poly_N = new int[face_N];
}

void PolygonalMesh::setPolygonCount(int index, int n)
{
	poly_N[index] = n;
	face[index] = new int[n];
	face[index][0] = -1;
}

void PolygonalMesh::computeFaceNormal()
{
	if(normal_f == NULL)
		normal_f = new float[face_N][3];
	for(int i=0; i<face_N; i++){
		int n = poly_N[i];
		int *f = face[i];
		if(f[0] < 0)
			continue;
		double nor[3];
		nor[0] = nor[1] = nor[2] = 0;
		for(int j=0; j<n; j++){
			int i1 = f[j];
			int i2 = f[(j+1)%n];
			nor[0] += (vertex[i1][1] - vertex[i2][1])*(vertex[i1][2] + vertex[i2][2]);
			nor[1] += (vertex[i1][2] - vertex[i2][2])*(vertex[i1][0] + vertex[i2][0]);
			nor[2] += (vertex[i1][0] - vertex[i2][0])*(vertex[i1][1] + vertex[i2][1]);
		}
		double len = LENGTH(nor);
		if((float)len != 0){
			normal_f[i][0] = (float)(nor[0]/len);
			normal_f[i][1] = (float)(nor[1]/len);
			normal_f[i][2] = (float)(nor[2]/len);
		}
		else{
			normal_f[i][0] = normal_f[i][1] = normal_f[i][2] = 0;
		}
	}
}

void PolygonalMesh::computeNormal()
{
	if(normal == NULL)
		normal = new float[vertex_N][3];
	int i;
	for(i=0; i<vertex_N; i++)
		normal[i][0] = normal[i][1] = normal[i][2] = 0;

	for(i=0; i<face_N; i++){
		int n = poly_N[i];
		int *f = face[i];
		double nor[3];
		nor[0] = nor[1] = nor[2] = 0;
		if(f[0] < 0)
			continue;
		int j;
		for(j=0; j<n; j++){
			int i1 = f[j];
			int i2 = f[(j+1)%n];
			nor[0] += (vertex[i1][1] - vertex[i2][1])*(vertex[i1][2] + vertex[i2][2]);
			nor[1] += (vertex[i1][2] - vertex[i2][2])*(vertex[i1][0] + vertex[i2][0]);
			nor[2] += (vertex[i1][0] - vertex[i2][0])*(vertex[i1][1] + vertex[i2][1]);
		}
		for(j=0; j<n; j++){
			int v = face[i][j];
			normal[v][0] += (float)nor[0];
			normal[v][1] += (float)nor[1];
			normal[v][2] += (float)nor[2];
		}
	}

	for(i=0; i<vertex_N; i++){
		double len = LENGTH(normal[i]);
		if((float)len != 0){
			normal[i][0] /= (float)len;
			normal[i][1] /= (float)len;
			normal[i][2] /= (float)len;
		}
	}
}

void PolygonalMesh::computeCenter()
{
	if(center == NULL)
		center = new float[face_N][3];
	for(int i=0; i<face_N; i++){
		center[i][0] = center[i][1] = center[i][2] = 0;
		int m = poly_N[i];
		for(int j=0; j<m; j++){
			center[i][0] += vertex[face[i][j]][0];
			center[i][1] += vertex[face[i][j]][1];
			center[i][2] += vertex[face[i][j]][2];
		}
		center[i][0] /= m;
		center[i][1] /= m;
		center[i][2] /= m;
	}
}


void PolygonalMesh::getBound(float &xmin, float &xmax, float &ymin, 
							 float &ymax, float &zmin, float &zmax)
{
	xmin = vertex[0][0];
	xmax = vertex[0][0];
	ymin = vertex[0][1];
	ymax = vertex[0][1];
	zmin = vertex[0][2];
	zmax = vertex[0][2];

	for(int i=1; i<vertex_N; i++){
		if(vertex[i][0] > xmax)
			xmax = vertex[i][0];
		if(vertex[i][0] < xmin)
			xmin = vertex[i][0];

		if(vertex[i][1] > ymax)
			ymax = vertex[i][1];
		if(vertex[i][1] < ymin)
			ymin = vertex[i][1];

		if(vertex[i][2] > zmax)
			zmax = vertex[i][2];
		if(vertex[i][2] < zmin)
			zmin = vertex[i][2];
	}	
}

float PolygonalMesh::rescale(float ori[3], float scale)
{
	float max[3], min[3];
	getBound(min[0], max[0], min[1], max[1], min[2], max[2]);

	float mx = 0.5f*(max[0]+min[0]);
	float my = 0.5f*(max[1]+min[1]);
	float mz = 0.5f*(max[2]+min[2]);

	float sx = max[0] - min[0];
	float sy = max[1] - min[1];
	float sz = max[2] - min[2];
	float s = scale/(float)sqrt(sx*sx + sy*sy + sz*sz);

	int i;
	for(i=0; i<vertex_N; i++){
	  float* p = vertex[i];
	  p[0] = (p[0]-mx)*s;
	  p[1] = (p[1]-my)*s;
	  p[2] = (p[2]-mz)*s;
	}

	ori[0] = mx;	ori[1] = my;	ori[2] = mz;
	return s;
}

float PolygonalMesh::fitIntoBox(float ct[3], float boxSize)
{
	float xMin, xMax, yMin, yMax, zMin, zMax;
	getBound(xMin, xMax, yMin, yMax, zMin, zMax);
//	printf("Original bounding box: (%f, %f, %f)-(%f, %f, %f)\n", xMin, yMin, zMin, xMax, yMax, zMax);
	//  calculate scale to fit into a 2*2*2 box
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
	int numVert = vertex_N;
	for (int i = 0; i < numVert; i++){
	//  pointer to vertex information
	float *p = vertex[i];
	p[0] = (p[0]-cX)*scale;
	p[1] = (p[1]-cY)*scale;
	p[2] = (p[2]-cZ)*scale;
	}

	return scale;
}


void PolygonalMesh::distanceFromPts(double &maxDis, double &avgDis, PointSet *pts)
{
	PQP_Model *pqp_meshModel = new PQP_Model();
	pqp_meshModel->BeginModel();
	int i;
	PQP_REAL p1[3], p2[3], p3[3];
	for(i = 0; i< face_N ; i++)
	{
		int vid = face[i][0];
		p1[0] = (PQP_REAL)(vertex[vid][0]);
		p1[1] = (PQP_REAL)(vertex[vid][1]);
		p1[2] = (PQP_REAL)(vertex[vid][2]);
		vid = face[i][1];
		p2[0] = (PQP_REAL)(vertex[vid][0]);
		p2[1] = (PQP_REAL)(vertex[vid][1]);
		p2[2] = (PQP_REAL)(vertex[vid][2]);
		vid = face[i][2];
		p3[0] = (PQP_REAL)(vertex[vid][0]);
		p3[1] = (PQP_REAL)(vertex[vid][1]);
		p3[2] = (PQP_REAL)(vertex[vid][2]);
		pqp_meshModel->AddTri(p1,p2,p3,i);
	}
	pqp_meshModel->EndModel();

	PQP_DistanceResult dres;	dres.last_tri = pqp_meshModel->last_tri;
	double maxD = 0.0;
	double avgD = 0.0;
	for(i = 0; i< pts->point_N; i++)
	{
	//	find minimal distance from a point to the mesh
		float *pt = pts->point[i];
		float *pn = pts->normal[i];
 		PQP_REAL p[3];
		p[0] = pt[0];	p[1] = pt[1];	p[2] = pt[2];
		PQP_Distance(&dres,pqp_meshModel, p,0.0,0.0);	
		
		if(maxD < dres.distance)
			maxD = dres.distance;
		avgD += dres.distance;
	}
	avgD = avgD/(double)pts->point_N;

	maxDis = maxD;	avgDis = avgD;

	delete pqp_meshModel;
}