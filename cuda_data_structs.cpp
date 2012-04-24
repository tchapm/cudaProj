/*
 *  cuda_data_struct.cpp
 *  cudaProj
 *
 *  Created by Tyler Chapman on 8/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include "cuda_data_structs.cuh"

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <fcntl.h>
#include <errno.h>
#include <math.h>
#include <sstream>

using namespace std;

float *cudaVector(int nelem)
{
	float *x=NULL;	
	if (nelem > 0)
		cudaMalloc((void**)&x, nelem*sizeof(float));
	return x;
}

float **cudaMatrix(int nx, int ny)
{
	float **device2DArray;
	float *h_temp[nx];
	cudaMalloc((void **)&device2DArray, nx*sizeof(float *));
	for (int i=0; i<nx; i++) {
		cudaMalloc((void **) &h_temp[i], ny*sizeof(float));
	}
	cudaMemcpy(device2DArray, h_temp, nx*sizeof(float *), cudaMemcpyHostToDevice);
	return device2DArray;
}

float ***cudaTensor(int nx, int ny, int nz)
{
	float ***x=NULL;
	long i, j;
	
	if (nx > 0 && ny > 0) 
		cudaMalloc((void**)&x, nx*sizeof(float **));
	for (i = 0; i < nx; ++i) {
		cudaMalloc((void**)&x[i], ny*sizeof(float *));
	}
	
	if (x == NULL){ 
		printf("Could not allocate float matrix with %d rows",nx);
		perror("errno");
		exit(1);
	}
	cudaMalloc((void**)&x[0][0], nx*ny*nz*sizeof(float));
	if (x[0] == NULL){ 
		printf("Could not allocate %d%d float matrix",nx,ny);
		perror("errno");
		exit(1);
		
	}
	for (i = 0; i < nx; ++i){
		for (j=0; j<ny; j++) {
			x[i][j] = x[0][0]+i*ny*nz+j*nz;
		}
	}
	return x;
}




