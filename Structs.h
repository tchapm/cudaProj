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
    int eGridSize;
	int tGridSize;
	int mGridSize;
    double* tempAxis;
    double* metalAxis;
    float depth; // in Mpc
	int nPixX;
    int nPixY;
	float pixToMpc;
	long binCenterSize;
    double theta;
    double phi;
    double epsilon;
    double a_ell;
    double b_ell;
    int n_ell;
    double rMax;
    float*** rebinnedCooling;
    float* metalArr;
    float* tempArr;
    float* emmArr;
	float accuracy;
	int nx;
	int ny;
	int nz;
	
};

#endif
