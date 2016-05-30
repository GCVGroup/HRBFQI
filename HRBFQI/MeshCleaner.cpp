#include "StdAfx.h"
#include "MeshCleaner.h"

CMeshCleaner::CMeshCleaner(void)
{
}

CMeshCleaner::~CMeshCleaner(void)
{
}

PolygonalMesh * CMeshCleaner::removeLowConfidenceGeometry(PolygonalMesh *oriMesh, HRBF *func, float confThreshold)
{
	if(oriMesh == NULL || func == NULL)
		return NULL;

	float thresholdConf = confThreshold;

	int vertex_N = oriMesh->vertex_N;
	float *conf = new float[vertex_N];
	memset(conf, 0, sizeof(float)*vertex_N);
	int i;

	//	computer the confidence for each vertex
	//-------------------------------------------
	//	only for rbfs with fixed radii

	for(i = 0; i< vertex_N; i++)
	{
		float c[3];
		c[0] = oriMesh->vertex[i][0];	c[1] = oriMesh->vertex[i][1]; c[2] = oriMesh->vertex[i][2];
		int size;
		int* t;
		t = func->tree->getIndexTable(size, c[0], c[1], c[2], func->support);
		int j;
		for(j = 0; j< size; j++)
		{
			int id = t[j];

			double r2 =( func->ps->point[id][0]-c[0])*( func->ps->point[id][0]-c[0])+( func->ps->point[id][1]-c[1])*(func->ps->point[id][1]-c[1])
				+( func->ps->point[id][2]-c[2])*(func->ps->point[id][2]-c[2]);
			conf[i] += func->weight(r2);
		}		
	}

	//---------------------------------------------------

	//	determine the invalid vertices and faces
	int face_N = oriMesh->face_N;
	bool *invalid_faces = new bool[face_N];
	memset(invalid_faces, 0, sizeof(bool)*face_N);

	int valid_faceN = 0;
	for(i = 0; i< face_N; i++)
	{
		if(conf[oriMesh->face[i][0]]<thresholdConf || conf[oriMesh->face[i][1]]<thresholdConf || conf[oriMesh->face[i][2]]<thresholdConf)
		{
			invalid_faces[i] = true;
		}
		else
			valid_faceN++;
	}

	int *id_corresponding = new int[vertex_N];
	int *move_ahead_number = new int[vertex_N];
	memset(move_ahead_number, 0, sizeof(int)*vertex_N);
	int move_number = 0;
	for(i = 0; i< vertex_N; i++)
	{
		move_ahead_number[i] = move_number;
		id_corresponding[i] = i-move_ahead_number[i];		
		if(conf[i]<thresholdConf)
		{
			move_number ++;
		}
	}
	delete[] move_ahead_number;
	
	//	update the mesh
	PolygonalMesh *newMesh = new PolygonalMesh();
	newMesh->setVertexCount(vertex_N-move_number);
	float (*vertex)[3] = newMesh->vertex;

	newMesh->setFaceCount(valid_faceN);
	int *poly_N = newMesh->poly_N;
	int **face = newMesh->face;

	int count = 0;
	for(i = 0; i< vertex_N; i++)
	{
		if(conf[i]<thresholdConf)
			continue;

		vertex[count][0] = oriMesh->vertex[i][0];
		vertex[count][1] = oriMesh->vertex[i][1];
		vertex[count][2] = oriMesh->vertex[i][2];

		count++;
	}

	count = 0;
	for(i = 0; i< face_N; i++)
	{
		if(!invalid_faces[i])
		{
			newMesh->setPolygonCount(count, 3);
			int *f = face[count];
			f[0] = id_corresponding[oriMesh->face[i][0]];
			f[1] = id_corresponding[oriMesh->face[i][1]];
			f[2] = id_corresponding[oriMesh->face[i][2]];

			count++;
		}
	}

	delete[] conf;
	delete[] invalid_faces;
	delete[] id_corresponding;

	return newMesh;
}

PolygonalMesh * CMeshCleaner::removeSmallIsolatedComponent(PolygonalMesh *oriMesh, int componentSize)
{
	oriMesh->setupTopologyInfo();

	vector<INTVECTOR> patches;

	int i;
	
	int vertexN = oriMesh->vertex_N;
	int faceN = oriMesh->face_N;

	bool *isVisited = new bool[vertexN];
	memset(isVisited, 0, vertexN*sizeof(bool));


	INTVECTOR seeds;
	for(i = 0; i< vertexN; i++)
	{
		if(isVisited[i])
			continue;

		seeds.clear();
		seeds.push_back(i);
		isVisited[i] = true;

		INTVECTOR patch;
		patch.push_back(i);

		while(!seeds.empty())
		{
			int seedId = seeds.back();
			seeds.pop_back();


			INTVECTOR edgeIds = oriMesh->m_vertLinkInfo[seedId].m_linkEdgeIds;

			int j;
			int edgeNum = edgeIds.size();
			for(j = 0; j< edgeNum; j++)
			{
				int eid = edgeIds[j];
				int vid1 = oriMesh->m_edgeLinkInfo[eid].m_linkVertIds[0];
				int vid2 = oriMesh->m_edgeLinkInfo[eid].m_linkVertIds[1];

				int vid = (vid1 == seedId ? vid2: vid1);
				if(!isVisited[vid])
				{
					seeds.push_back(vid);
					patch.push_back(vid);
					isVisited[vid] = true;
				}
			}
		}

		patches.push_back(patch);
	}
	
	//---------------------------------------------------

	//	determine the invalid vertices and faces
	bool *invalid_vertices;
	invalid_vertices = isVisited;
	memset(invalid_vertices, 0, vertexN*sizeof(bool));


	int patchNum = patches.size();
//	cout<<"Patches: "<< patchNum<< endl;
	for(i = 0; i< patchNum; i++)
	{
		INTVECTOR patch;
		patch = patches[i];
		int patchSize = patch.size();
//		cout<<i<<"patch size: "<< patchSize<< endl;
		if(patchSize < componentSize)
		{
			int j; 
			for(j = 0; j< patchSize; j++)
			{
				int vid = patch[j];
				invalid_vertices[vid] = true;
			}
		}
	}


	bool *invalid_faces = new bool[faceN];
	memset(invalid_faces, 0, sizeof(bool)*faceN);

	int valid_faceN = 0;
	for(i = 0; i< faceN; i++)
	{
		if(invalid_vertices[oriMesh->face[i][0]] || invalid_vertices[oriMesh->face[i][1]] || invalid_vertices[oriMesh->face[i][2]])
		{
			invalid_faces[i] = true;
		}
		else
			valid_faceN++;
	}

	int *id_corresponding = new int[vertexN];
	int *move_ahead_number = new int[vertexN];
	memset(move_ahead_number, 0, sizeof(int)*vertexN);
	int move_number = 0;
	for(i = 0; i< vertexN; i++)
	{
		move_ahead_number[i] = move_number;
		id_corresponding[i] = i-move_ahead_number[i];		
		if(invalid_vertices[i])
		{
			move_number ++;
		}
	}

	delete[] move_ahead_number;
	
	//	update the mesh
	PolygonalMesh *newMesh = new PolygonalMesh();
	newMesh->setVertexCount(vertexN-move_number);
	float (*vertex)[3] = newMesh->vertex;

	newMesh->setFaceCount(valid_faceN);
	int *poly_N = newMesh->poly_N;
	int **face = newMesh->face;

	int count = 0;
	for(i = 0; i< vertexN; i++)
	{
		if(invalid_vertices[i])
			continue;

		vertex[count][0] = oriMesh->vertex[i][0];
		vertex[count][1] = oriMesh->vertex[i][1];
		vertex[count][2] = oriMesh->vertex[i][2];

		count++;
	}

	count = 0;
	for(i = 0; i< faceN; i++)
	{
		if(invalid_faces[i])
			continue;

		newMesh->setPolygonCount(count, 3);
		int *f = face[count];
		f[0] = id_corresponding[oriMesh->face[i][0]];
		f[1] = id_corresponding[oriMesh->face[i][1]];
		f[2] = id_corresponding[oriMesh->face[i][2]];

		count++;
	}

	delete[] invalid_vertices;
	delete[] invalid_faces;
	delete[] id_corresponding;


	return newMesh;
}
