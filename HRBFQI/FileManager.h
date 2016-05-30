// FileManager.h: FileManager 
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FILEMANAGER_H__84206D95_38F8_41BC_BA6F_F12B3D1109F2__INCLUDED_)
#define AFX_FILEMANAGER_H__84206D95_38F8_41BC_BA6F_F12B3D1109F2__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "../DataStructure/PolygonalMesh.h"
#include "../DataStructure/PointSet.h"
#include "stdio.h"

class FileManager  
{
protected:
	void saveWaveFront(PolygonalMesh *mesh, float oriCt[3], float scale = 1.0f);

	FILE* file;
	CString file_name;
	CString file_ext;

public:
	void save(PolygonalMesh *mesh, float oriCt[3], float scale = 1.0f);
	void setFile(FILE* file, CString file_name, CString file_ext);
	void open(PointSet* ps);
	FileManager();
	virtual ~FileManager();

};

#endif // !defined(AFX_FILEMANAGER_H__84206D95_38F8_41BC_BA6F_F12B3D1109F2__INCLUDED_)
