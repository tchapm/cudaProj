//
//  Structs.h
//  cudaProj
//
//  Created by Tyler Chapman on 3/6/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef cudaProj_Structs_h
#define cudaProj_Structs_h

struct constants {
	float depth; // in Mpc
	int nPixX;
    int nPixY;
	float pixToMpc;
	long binCenterSize;
	float accuracy;
	int nGrid;
	int nx;
	int ny;
	int nz;
	int eGridSize;
	int tGridSize;
	int mGridSize;
    double *tempAxis;
    double *metalAxis;
	
};

#endif
