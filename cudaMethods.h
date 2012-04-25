//
//  cudaMethods.h
//  cudaProj
//
//  Created by Tyler Chapman on 2/21/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef cudaProj_cudaMethods_h
#define cudaProj_cudaMethods_h

__global__ void cudFunctionTest(double *testArr_d);

__device__ void modifyTensor(double* oneDArray, int index);

__device__ double tenRetrieveD(double* oneDArray, int nx, int ny, int nz, int x, int y, int z);

__device__ double linearInterpZ(double* interpMat, int x, int y, double z, constants theConst);

__device__ double linearInterpEll(double* interpMat, double ell, constants theConst);

__device__ int getLowerIndex(double* axisArr, double value, int arrLength);

__device__ double bilinInterpVal(double* interpMat, double y, double z, int energyBin, double* tempGrid, double* metalGrid, constants theConst);


__device__ double getEll(int x, int y, int z, double* rotMat, double a, double b);

__global__ void integrate2(double* rebinCool, double* tempArr, double* metalArr, double* emmArr, double* integral, double* tempGrid, double* metalGrid, double* rotMat, constants theConst, bool debugging);

__global__ void integrate(double* rebinCool, double* tempArr, double* metalArr, double* emmArr, double* integral, double* tempGrid, double* metalGrid, constants theConst, bool debugging);

__global__ void tempInit(void(*f)(double), double * temperature, constants theConst);

__device__ double function (double rNot, double l) ;

__device__ void function2 (double x);
#endif
