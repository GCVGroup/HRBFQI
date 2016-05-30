// HRBFQI.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <string>
#include <time.h>
#include "../DataStructure/PointSet.h"
#include "../DataStructure/PolygonalMesh.h"
#include "FileManager.h"
#include "../polygonizer/Polygonizer.h"
#include "../DataStructure/OctTree.h"
#include "HRBF.h"
#include "MeshCleaner.h"
#include "omp.h"

#define PI 3.14159265

using namespace std;

int main(int argc, char* argv[])
{
    std::cout << "Hermite RBF-based Reconstruction "
        << "(built on " __DATE__ ", " __TIME__ ")." << std::endl;


	if(argc != 10)
	{
		cout<<"HBRFQI inputfile.pwn(.pwnc) outputfile.obj <if_rescale> <if_output_difference_gradient_normal> <support_size_scale>"<<
			" <parameter_eta> <isosurface_extract_stepsize> <confidence_threshold> <component_size>"<<endl;
		return 1;
	}	

	//	set parameter values
	double supportSizeScale, etaValue, gridEdgeLen, confThreshold;
	bool ifRescale, ifOutputDifferenceGradientNormal;
	int componentSize;
	
	ifRescale = atoi(argv[3]);
	ifOutputDifferenceGradientNormal = atoi(argv[4]);
	supportSizeScale = atof(argv[5]);
	etaValue = atof(argv[6]);
	gridEdgeLen = atof(argv[7]);
	confThreshold = atof(argv[8]);
	componentSize = atoi(argv[9]);

	//-------------------------------------------------------------------------
	//	read file
	FileManager * file_manager = new FileManager;
	PointSet *ps = new PointSet;
	float fScale = 1.0f;
	float fCenter[3] = {0.0f, 0.0f, 0.0f};
	cout<<"********************************************************************"<<endl;
	cout<<"Read points file..."<<argv[1]<<endl;
	char input_filename[256];
	strcpy(input_filename,"Data\\");
	strcat(input_filename,argv[1]);
//	int length=strlen(filename);
	FILE* file;
	if(strstr(argv[1], ".pwn") != NULL)
	{
		file = fopen(input_filename,"r");
		file_manager->setFile(file, input_filename, "pwn");
		file_manager->open(ps);
		fclose(file);
	}
	else
	{
		cout<<"Input file should be a pwn file."<<endl;
		exit(0);
	}
	if(ifRescale)
	{
		fScale = ps->fitIntoBox(fCenter, 1.0f);
	}
	cout<<argv[1]<<": file input finished!"<<endl;
	cout<<ps->point_N<<" points input"<<endl;


	//-------------------------------------------------------------------------
	//	Interpolation or Quasi-interpolation
	cout<<endl<<"********************************************************************"<<endl;
	cout<<"Fitting..."<<endl;
	clock_t start = clock();

	HRBF *hrbf = new HRBF();
	hrbf->setPointSet(ps);

	//Estimation of recommended support size
	float T = 0.75*hrbf->getAveragedLeafSize();
//	cout<<"The initial support radius: "<< T <<endl;
	
	hrbf->support = supportSizeScale*T;
	T = hrbf->support;

	float normal_smooth;

	normal_smooth = etaValue;///(T*T);//100.0f/(T*T);
	hrbf->fit(hrbf->support, normal_smooth);
	cout<<"Done!"<<endl;
	if(!ifOutputDifferenceGradientNormal)
	{
	cout<<"The time for fitting(including the octree construction): "<< (float)(clock()-start)/ CLOCKS_PER_SEC<<"s;"<<endl;
	//---------------
	//	for debug
	//cout<<"--------------------------"<<endl;
	//cout<<"The final support radius: "<< hrbf->support<<endl;
	//int neighborNum = hrbf->getMaximalNeighborsInSupport(hrbf->support);	//	Compute a maximal neighbor number according to the given support size
	//cout<<"The maximal neighbor number: "<< neighborNum<< endl;
	//float rho_min = (5.0*neighborNum+sqrt(25.0*neighborNum*neighborNum+2240.0*(1+normal_smooth)))/(8.0*(1.0+normal_smooth)); 
	//cout<<"The parameter: "<< normal_smooth<< endl;//fit_dia->m_normal_smooth
	//cout<<"The minimal support radius: "<< rho_min<<endl;
	}
	//-------------------------------------------------------------------------
	//	Iso-surface extraction
	cout<<endl<<"********************************************************************"<<endl;
	cout<<"Surface extraction..."<<endl;
	start = clock();
	float x_max, x_min, y_max, y_min, z_max, z_min;
	ps->getBound(x_min, x_max, y_min, y_max, z_min, z_max);
	float bound_X = x_max - x_min;
	float bound_Y = y_max - y_min;
	float bound_Z = z_max - z_min;
	float grid_size = gridEdgeLen; 

	float space = grid_size;

	Polygonizer* poly = new Polygonizer;

	poly->spaceX = space;
	poly->spaceY = space;
	poly->spaceZ = space;

	poly->originX = x_min - 5*space;
	poly->originY = y_min - 5*space;
	poly->originZ = z_min - 5*space;

	poly->dimX = (int)((bound_X)/space) + 10;
	poly->dimY = (int)((bound_Y)/space) + 10;
	poly->dimZ = (int)((bound_Z)/space) + 10;

	poly->func = hrbf;			

	PolygonalMesh *mesh;
	mesh = NULL;
	mesh = poly->dualContouring(0.001f, 0.01f); 

	cout<<"Done!"<<endl;
	if(!ifOutputDifferenceGradientNormal)
	{
		cout<<"The time for getting triangles:"<<((float)(clock()-start))/ CLOCKS_PER_SEC<<"s;"<<endl;
		cout<<"The extracted mesh has " << mesh->vertex_N<<" vertices and "<< mesh->face_N<<" triangles."<<endl;
	}

	//---------------
	//	show the differences between gradients and normals
	if(ifOutputDifferenceGradientNormal)
	{
		double maxAngle, aveAngle;

		hrbf->differenceFuncGradientAndNormalPt(maxAngle, aveAngle);
		cout<<"The maximal angle between gradients and normals: "<<maxAngle<<endl;
		cout<<"The average angle between gradients and normals: "<< aveAngle << endl<< endl;

		// distances between points and the reconstructed mesh
		double maxDist, aveDist;
		mesh->distanceFromPts(maxDist, aveDist, ps);
		cout<<"The maximal distance between points and the reconstructed mesh: "<< maxDist<<endl;
		cout<<"The maximal distance between points and the reconstructed mesh: "<< aveDist << endl;

	}
	//-------------------------------------------------------------------------
	//	Mesh cleaning
	if(!ifOutputDifferenceGradientNormal)
	{

	if(confThreshold>0.0 || componentSize>0.0)
	{
		CMeshCleaner *meshCleaner = new CMeshCleaner;
		PolygonalMesh *newMesh;
		newMesh = NULL;

		cout<<endl<<"********************************************************************"<<endl;
		cout<<"Mesh cleaning..."<<endl;
		cout<<"--------------------------"<<endl;
		cout<<"Removing low-confidence geometry (threshold "<< confThreshold << ")..." <<endl;

		clock_t duration = 0;
		start = clock();
		newMesh = meshCleaner->removeLowConfidenceGeometry(mesh, hrbf, confThreshold);
		duration += clock()-start;

		cout<<"Removed "<< mesh->vertex_N - newMesh->vertex_N << " vertices and " << mesh->face_N - newMesh->face_N<< " triangles."<<endl;


		delete mesh;
		//delete ps;
		//delete hrbf;
		//delete poly;

		if (componentSize > 0)
		{
			cout<<"--------------------------"<<endl;
			cout<<"Removing isolated components with < " << componentSize << " vertices..." <<endl;
			
			mesh = newMesh;
			start = clock();
			newMesh = meshCleaner->removeSmallIsolatedComponent(mesh, componentSize);
			duration += clock()-start;

			cout<<"Removed "<< mesh->vertex_N - newMesh->vertex_N << " vertices and " << mesh->face_N - newMesh->face_N<< " triangles."<<endl;

			delete mesh;
		}
		cout<<"--------------------------"<<endl;
		cout<<"The time for mesh cleaning:"<<((float)(duration))/ CLOCKS_PER_SEC<<"s."<<endl;
		mesh = newMesh;
		delete meshCleaner;
	}
	//-------------------------------------------------------------------------
	//	write file
	cout<<endl<<"********************************************************************"<<endl;
	cout<<"Write mesh file..."<<endl;

	if(mesh == NULL || mesh->face_N == 0)
		return 1;
	char output_filename[256];
	strcpy(output_filename,"Data\\");
	strcat(output_filename,argv[2]);

	file = fopen(output_filename,"w");
	file_manager->setFile(file, output_filename, "obj");
//	float tempCenter[3] = {0,0,0};
	file_manager->save(mesh, fCenter, fScale);//mesh, tempCenter);

	fclose(file);
	cout<<argv[2]<<": file output finished!"<<endl;
	cout<<mesh->vertex_N <<" vertices and " << mesh->face_N << " triangles output."<<endl;
	cout<<endl<<"********************************************************************"<<endl;
	}

	delete mesh;
	delete ps;
	delete hrbf;
	delete poly;
	delete file_manager;
	
	cout<<"Press Enter to return!"<<endl;
	getchar();
    return 0;
}

