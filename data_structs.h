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
#include "Structs.h"
float **sci_fmatrix(int nx, int ny);

double **sci_fmatrixD(int nx, int ny);

float ***sci_ftensor(int nx, int ny, int nz);

double ***sci_ftensorD(int nx, int ny, int nz);

float *tensorTo1DArray(float ***tensor3D, int nx, int ny, int nz);

double *tensorTo1DArray(double ***tensor3D, int nx, int ny, int nz);

float *oneDtensorToMatrix(float *tensor1D, int nx, int ny, int nz, int z);

float *matrixTo1DArray(float **matrix2D, int nx, int ny);

double *metalArrInit(int nx, int ny, int nz);

double *metalArrInit(int nEll, double a_ell, double b_ell, double rMax);

double *emmArrInit(int nx, int ny, int nz);

double *emmArrInit(int nEll, double a_ell, double b_ell, double rMax);

double *tempArrInit(int nx, int ny, int nz);

double *tempArrInit(int nEll, double a_ell, double b_ell, double rMax);

double *energyArrInit(int rebinSize);

double *ellArrInit(int nEll, double rMax);

double *tempGridInit(int tGridSize, double *tempAxis);

double *metalGridInit(int mGridSize, double *metalAxis);

double *rotMatInit(double theta, double phi, double epsilon);

double ***makeIntegralMatrix(double *integral_h, int nPixX, int nPixY, int binCenterSize);

double oneDArrayToMatrix(double oneDArray[], int nx, int ny);
    
double get3DValFrom1DArray(double oneDArray[], int x, int y, int z, int nx, int ny, int nz);

double get2DValFrom1DArray(double oneDArray[], int x, int y, int nx, int ny);

double tenRetrieveH(double oneDArray[], int nx, int ny, int nz, int x, int y, int z);

void printSpectra(double*** integralMatrix, double* energyArr, constants theConst);

void sumSpectra(double*** integralMatrix, double* energyArr, constants theConst);

#endif
