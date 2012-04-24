//
//  cudaMethods.h
//  cudaProj
//
//  Created by Tyler Chapman on 2/21/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef cudaProj_cudaMethods_h
#define cudaProj_cudaMethods_h

__global__ void integrate(float* rebinCool, float* tempArr, float* metalArr, float* emmArr, float* integral, float* tempGrid, float* metalGrid, constants theConst, bool debugging);

__global__ void cudFunctionTest(float *testArr_d);

__device__ void modifyTensor(float* oneDArray, int index);

__device__ float tenRetrieveD(float* oneDArray, int nx, int ny, int nz, int x, int y, int z);

__device__ float linearInterpZ(float* interpMat, int x, int y, float z, constants theConst);

__device__ int getLowerIndex(float* axisArr, float value, int arrLength);

__device__ float bilinInterpVal(float* interpMat, float y, float z, int energyBin, float* tempGrid, float* metalGrid, constants theConst);

__global__ void tempInit(void(*f)(float), float * temperature, constants theConst);

__device__ float function (float rNot, float l) ;

__device__ void function2 (float x);
#endif
