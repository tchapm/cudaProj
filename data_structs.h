/*
 *  data_structs.h
 *  cudaProj
 *
 *  Created by Tyler Chapman on 8/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef DATA_STRUCTS_H
#define DATA_STRUCTS_H

float **sci_fmatrix(int nx, int ny);

float ***sci_ftensor(int nx, int ny, int nz);

float *tensorTo1DArray(float ***tensor3D, int nx, int ny, int nz);

float *matrixTo1DArray(float **matrix2D, int nx, int ny);

float *metalArrInit(int nx, int ny, int nz);

float *emmArrInit(int nx, int ny, int nz);

float *tempArrInit(int nx, int ny, int nz);

float *energyArrInit(int rebinSize);

float *tempGridInit(int tGridSize, double *tempAxis);

float *metalGridInit(int mGridSize, double *metalAxis);

float ***makeIntegralMatrix(float *integral_h, int nPixX, int nPixY, int binCenterSize);

float get3DValFrom1DArray(float oneDArray[], int x, int y, int z, int nx, int ny, int nz);

float get2DValFrom1DArray(float oneDArray[], int x, int y, int nx, int ny);

float tenRetrieveH(float oneDArray[], int nx, int ny, int nz, int x, int y, int z);

void printSpectra(float *** integralMatrix, int nPixX, int nPixY);

#endif
