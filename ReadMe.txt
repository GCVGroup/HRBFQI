The source code is for the method presented in

Shengjun Liu, Charlie C.L. Wang, Guido Brunnett, Jun Wang. "A Closed-Form Formulation of HRBF-Based Surface Reconstruction by Approximate Solution". Computer-aided Design, accepted, 2016. (Special issue for SPM2016).


It can be compiled with Microsoft Visual Studio 2015, and run on the operating system --- Windows 10.





*****************************************************************
  	    Hermite RBF Quasi-interpolation (HRBFQI)
*****************************************************************
By: Shengjun Liu
email: shjliu.cg@csu.edu.cn
webpage: github.com/GCVGroup/HRBFQI
*****************************************************************

1.Copyright
-----------

- HRBFQI is developed by Shengjun Liu for research use. All rights about the program (esp. surface reconstruction) are reserved by Shengjun Liu. This C++ source codes are available only to a primary user for academic purposes. No secondary use, such as copy, distribution, diversion, business purpose, etc., is allowed. In no event shall the author be liable to any party for direct, indirect, special, incidental, or consequential damage arising out of the use of this program. HRBFQI is self-contained. 


2.Download
----------

- The source code (MSVS 2015), as well as the testing data, from the page: 
  
  https://github.com/GCVGroup/HRBFQI.

3.Installing & Compiling (Windows+MSVS2015)
-------------------------------------------

- Simply download the source code to a suitable place and use MSVS2015 to building the project. It can be compiled with Microsoft Visual Studio 2015, and run on the operating system --- Windows 10.


4.Usage
-------

- After the compilation you may try the tool HRBFQI which is inside the ./bin directory:

*Note* The data files should be located in ./bin/data directory.

- In console mode, the command will be:

   HBRFQI.ext inputfile.pwn outputfile.obj IR IODGN SSS PE IES CT CS

- Parameters:

  inputfile.pwn: the input file should be a pwn file.

  outputfile.obj: the output file should be a obj file.

  IR: (I)f_(R)escale. It is a bool parameter for determining if scaling the points into a <-1,-1,-1>*<1,1,1> box. 1 for scaling and 0 for not. 
  
  IODGN: (I)f_(O)utput_(D)ifference_(G)radient_(N)ormal. It is a bool value for showing the differences between reconstructed surface gradients and given normals at the points, and the distances from points to the reconstructed surface, or not. If IODGN is true, it means the program only shows the information, and the extracted mesh will not be output.

  SSS: (S)upport_(S)ize_(S)cale. SSS is a scalor to the initial support size computed with the given data.

  PE: (P)arameter_(E)ta. It is the parameter for regularization.

  IES: (I)sosurface_(E)xtract_(S)stepsize. IES is the size of the grid edge for isosurface extraction.

  CT: (C)onfidence_(T)hreshold. CT is a threshold for removing the triangles in which there is at least one vertex with potential value less than CT.

  CS: (C)omponent_(S)ize. It is a threshold for removing the patches with less than CS vertices.

- Examples:

  There are many examples provided in the ./bin directory. All data shown in the figures of our paper are provided. And .bat files in the ./bin directory can exacute directly to generate the results presented in figures or tables which are denoted as their file names.


5.File format
-------------

- pwn file: point with normal file (the files in the fold of \bin\data are with normal vectors inconsistently oriented).

#points (=N)
x1 y1 z1
x2 y2 z2
...
xN yN zN
nx1 ny1 nz1
nx2 ny2 nz2
...
nxN nyN nzN

- obj file: a general wavefront file with output mesh
