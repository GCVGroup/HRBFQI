// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#ifndef _WIN32_WINNT		// Allow use of features specific to Windows XP or later.                   
#define _WIN32_WINNT 0x0501	// Change this to the appropriate value to target other versions of Windows.
#endif						

//#include <stdio.h>
//#include <tchar.h>

#define VC_EXTRALEAN		// 

#include <afxwin.h>         // 
#include <afxext.h>         // 
#include <afxdisp.h>        // 
#include <afxdtctl.h>		// 
#ifndef _AFX_NO_AFXCMN_SUPPORT
#include <afxcmn.h>			// 
#endif // _AFX_NO_AFXCMN_SUPPORT

#define THREADS_NUM	8
#define MAX_NEIGHBOURS_NUM 10000

// TODO: reference additional headers your program requires here
