/*
 *  cuda_data_structs.h
 *  cudaProj
 *
 *  Created by Tyler Chapman on 8/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CUDA_DATA_STRUCTS_H
#define CUDA_DATA_STRUCTS_H

float *cudaVector(int nelem);

float **cudaMatrix(int nx, int ny);

float ***cudaTensor(int nx, int ny, int nz);

#endif