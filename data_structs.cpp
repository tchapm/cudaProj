/*
 *  data_structs.cpp
 *  cudaProj
 *
 *  Created by Tyler Chapman on 8/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "data_structs.h"

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <fcntl.h>
#include <errno.h>
#include <math.h>
#include <sstream>
//#include <sciutils.h>
//#include <jaco.h>

using namespace std;


float **sci_fmatrix(int nx, int ny)
{
	float **x=NULL;
	long i;
	
	if (nx > 0 && ny > 0) 
		x = (float **)malloc(nx*sizeof(float *));
	if (x == NULL){ 
		printf("Could not allocate float matrix with %d rows",nx);
		perror("errno");
		exit(1);
	}
	
	x[0] = (float *)malloc(nx*ny*sizeof(float));
	if (x[0] == NULL){ 
		printf("Could not allocate %d%d float matrix",nx,ny);
		perror("errno");
		exit(1);
		
	}
	for (i = 1; i < nx; ++i) 
		x[i] = x[0]+i*ny;
	
	return x;
}

double ***sci_ftensorD(int nx, int ny, int nz)
{
	double ***x=NULL;
	long i, j;
	
	if (nx > 0 && ny > 0) 
		x = (double ***)malloc(nx*sizeof(double **));
	
	for (i = 0; i < nx; ++i) {
		x[i] = (double **)malloc(ny*sizeof(double *));
	}
	
	if (x == NULL){ 
		printf("Could not allocate float matrix with %d rows",nx);
		perror("errno");
		exit(1);
	}
	
	x[0][0] = (double *)malloc(nx*ny*nz*sizeof(double));
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
float ***sci_ftensor(int nx, int ny, int nz)
{
	float ***x=NULL;
	long i, j;
	
	if (nx > 0 && ny > 0) 
		x = (float ***)malloc(nx*sizeof(float **));
	
	for (i = 0; i < nx; ++i) {
		x[i] = (float **)malloc(ny*sizeof(float *));
	}
	
	if (x == NULL){ 
		printf("Could not allocate float matrix with %d rows",nx);
		perror("errno");
		exit(1);
	}
	
	x[0][0] = (float *)malloc(nx*ny*nz*sizeof(float));
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

float *tensorTo1DArray(float ***tensor3D, int nx, int ny, int nz){
	float *oneDArray = new float[nx*ny*nz];
	
	int arrIndex = 0;
	for (int i=0; i<nx; i++) {
		for (int j=0; j<ny; j++) {
			for (int k=0; k<nz; k++) {
				oneDArray[arrIndex]=tensor3D[i][j][k];
				if (i==150 && j==2 && k==3) {
					//cout << "\noneDArr[150][2][3] = " << oneDArray[arrIndex] << "\tindex = "  << arrIndex << endl;
				}
				arrIndex++;
			}
		}
	}
	return oneDArray;
}

float *tensorTo1DArray(double ***tensor3D, int nx, int ny, int nz){
	float *oneDArray = new float[nx*ny*nz];
	
	int arrIndex = 0;
	for (int i=0; i<nx; i++) {
		for (int j=0; j<ny; j++) {
			for (int k=0; k<nz; k++) {
				oneDArray[arrIndex]=(float)tensor3D[i][j][k];
				if (i==150 && j==2 && k==3) {
					//cout << "\noneDArr[150][2][3] = " << oneDArray[arrIndex] << "\tindex = "  << arrIndex << endl;
				}
				arrIndex++;
			}
		}
	}
	return oneDArray;
}

float *oneDtensorToMatrix(float *tensor1D, int nx, int ny, int nz, int z){
	float *oneDMatrix = new float[nx*ny];
	
	for (int m=0; m<nx; m++) {
		for (int n=0; n<ny; n++) {
			oneDMatrix[m*ny + n] = tensor1D[m*ny*nz + n*nz + z];
		}
	}
	return oneDMatrix;
}

float *matrixTo1DArray(float **matrix2D, int nx, int ny){
	float *oneDArray = new float[nx*ny];
	int arrIndex = 0;
	for (int i=0; i<nx; i++) {
		for (int j=0; j<ny; j++) {
			oneDArray[arrIndex]=matrix2D[i][j];
			if (i==1 && j==2) {
				//cout << arrIndex;
			}
			arrIndex++;
		}
	}
	return oneDArray;
}

float *metalArrInit(int nx, int ny, int nz){
	int x, y, z;
	float r;
	float *metalArr = new float[nx*ny*nz];
    
	for (x=0; x<nx; x++) {
		for(y=0; y<ny; y++){
			for (z=0; z<nz; z++) {
				r=powf((powf(x-(nx/2-1),2)+powf(y-(ny/2-1),2)+powf(z-(nz/2-1), 2)), 0.5);
				metalArr[x*ny*nz + y*nz + z] = log10(1.0);//log10(powf((r+1), -0.5));
				if (x ==1 && y == 2 && z==3) {
					//printf("metalTen[1][2][3] = %f\t", powf((r+1), -0.5));
				}
			}
		}
	}
	return metalArr;
}

float *emmArrInit(int nx, int ny, int nz){
	int x, y, z;
	float r;
	float *emmArr = new float[nx*ny*nz];
    float l, A, B, C, rNot;
    A = 1.0;
    B = 1.0;
    C = 1.0;
    rNot = 0.5;
	for (x=0; x<nx; x++) {
		for(y=0; y<ny; y++){
			for (z=0; z<nz; z++) {
                l =powf((powf(float(x-(nx/2-1))/A,2)+powf(float(y-(ny/2-1))/B,2)+powf(float(z-(nz/2-1))/C,2)), 0.5);
				r=powf((powf(x-(nx/2-1),2)+powf(y-(ny/2-1),2)+powf(z-(nz/2-1), 2)), 0.5);
				emmArr[x*ny*nz + y*nz + z] = log10(powf(rNot+0.01*l,-3));//1.0;//powf((pow(r,2)+1), -3);
				if (x ==1 && y == 2 && z==3) {
					//printf("emmTen[1][2][3] = %f\t", powf((pow(r,2)+1), -3));
				}
			}
		}
	}
	return emmArr;
}

float *tempArrInit(int nx, int ny, int nz){
	int x, y, z;
	float r;
	float *tempArr = new float[nx*ny*nz];
    //    float max = 0;
    float l, A, B, C, rNot;
    A = 1.0;
    B = 1.0;
    C = 1.0;
    rNot = 0.5; //500 kPc * 0.01 pixel/Mpc
	for (x=0; x<nx; x++) {
		for(y=0; y<ny; y++){
			for (z=0; z<nz; z++) {
                //need to make l^2 makes and error when done
                l =powf((powf(float(x-(nx/2-1))/A,2)+powf(float(y-(ny/2-1))/B,2)+powf(float(z-(nz/2-1))/C,2)), 0.5);
				r=powf((powf(x-(nx/2-1),2)+powf(y-(ny/2-1),2)+powf(z-(nz/2-1), 2)), 0.5);
				tempArr[x*ny*nz + y*nz + z]=log10(6.7*powf(rNot+0.01*l,-1.0));//log10(2.0);//log10(6*powf((r+1), -0.5)); //Kev 11600000
                //                printf("tempTen[1][2][3] = %f\t",tempArr[x*ny*nz + y*nz + z]);log10((1/4)*
                //                if (x ==1 && y == 2 && z==3) {
                //                    printf("tempTen[1][2][3] = %f\t", tempArr[x*ny*nz + y*nz + z]);
                //                    printf("l = %f\n",l);
                //                }
                //				}
                //                if (tempArr[x*ny*nz + y*nz + z]>max){
                //                    max = tempArr[x*ny*nz + y*nz + z];
                //                    xmax = x;ymax = y;zmax =z;
                //                }
			}
		}
	}
    //    printf("\nmax tempInit = %f at [%d][%d][%d]\n ", max,xmax,ymax,zmax);
	return tempArr;
}

float *tempArrInit(int nEll, double a_ell, double b_ell, double rMax){
	float *tempArr = new float[nEll];
    //    float max = 0;
    float rNot = 0.5; //500 kPc * 0.01 pixel/Mpc
    float *ellArr = new float[nEll];
    ellArr = ellArrInit(nEll, rMax);
    for(int k=0; k<nEll; k++){
        tempArr[k]=log10(6.7*powf(rNot+0.01*ellArr[k],-1.0));
    }

	return tempArr;
}

float *emmArrInit(int nEll, double a_ell, double b_ell, double rMax){
	float *emmArr = new float[nEll];
    //    float max = 0;
    float rNot = 0.5; //500 kPc * 0.01 pixel/Mpc
    float *ellArr = new float[nEll];
    ellArr = ellArrInit(nEll, rMax);
    for(int k=0; k<nEll; k++){
        emmArr[k] = log10(powf(rNot+0.01*ellArr[k],-3));
    }
	return emmArr;
}

float *metalArrInit(int nEll, double a_ell, double b_ell, double rMax){
	float *metalArr = new float[nEll];
    float *ellArr = ellArrInit(nEll, rMax);
    for(int k=0; k<nEll; k++){
        metalArr[k] = log10(1.0);
    }
	return metalArr;
}

float *ellArrInit(int nEll, double rMax){
	int x;
	float eMin = 0.0; //in Kev
	float eMax = rMax;
	float binWidth = (eMax-eMin)/nEll;
	float *ellArr = new float[nEll];
	ellArr[0] = eMin;
	for (x=1; x<rMax; x++) {
		ellArr[x] = ellArr[x-1]+binWidth;
	}	
	return ellArr;
}
float *energyArrInit(int rebinSize){
	int x;
	float eMin = 0.2; //in Kev
	float eMax = 10.0;
	float binWidth = (eMax-eMin)/rebinSize;
	float *energyArr = new float[rebinSize];
	energyArr[0] = eMin;
	for (x=1; x<rebinSize; x++) {
		energyArr[x] = energyArr[x-1]+binWidth;
	}
	return energyArr;
	
}

float *tempGridInit(int tGridSize, double *tempAxis){
    float *tempGrid_h = (float*) calloc(tGridSize, sizeof(float));
    for(int i =0; i<tGridSize; i++){
		tempGrid_h[i] = pow(10,tempAxis[i]);
	}
    return tempGrid_h;
}

float *metalGridInit(int mGridSize, double *metalAxis){
    float *metalGrid_h = (float*) calloc(mGridSize, sizeof(float));
    for(int i =0; i<mGridSize; i++){
		metalGrid_h[i] = pow(10,metalAxis[i]);
	}
    return metalGrid_h;
}

float *rotMatInit(double theta, double phi, double epsilon){
    float *rotMat = new float[9];
    rotMat[0] = cos(epsilon)*cos(phi)-cos(theta)*sin(phi)*sin(epsilon);
    rotMat[1] = cos(epsilon)*sin(phi)+cos(theta)*cos(phi)*sin(epsilon);
    rotMat[2] = sin(epsilon)*sin(theta);
    rotMat[3] = -sin(epsilon)*cos(phi)-cos(theta)*sin(phi)*cos(epsilon);
    rotMat[4] = -sin(epsilon)*sin(phi)+cos(theta)*cos(phi)*cos(epsilon);
    rotMat[5] = cos(epsilon)*sin(theta);
    rotMat[6] = sin(theta)*sin(phi);
    rotMat[7] = -sin(theta)*cos(phi);
    rotMat[8] = cos(theta);
    return rotMat;
}

float ***makeIntegralMatrix(float *integral_h, int nPixX, int nPixY, int binCenterSize){
    float ***integralMatrix;
    integralMatrix = sci_ftensor(nPixX, nPixY, binCenterSize);
    for (int x=0; x<nPixX; x++) {
        for (int y=0; y<nPixY; y++) {
            for (int z=0; z<binCenterSize; z++) {
                integralMatrix[x][y][z] = log10(tenRetrieveH(integral_h, nPixX, nPixY, binCenterSize, x, y, z));
            }
        }
    }
    return integralMatrix;
}

float oneDArrayToMatrix(float oneDArray[], int nx, int ny){
    /*	int x, y, z;
     for (x=0; x<nx; x++) {
     for(y=0; y<ny; y++){
     for (z=0; z<nz; z++) { */
	return oneDArray[0];
    
}

float get3DValFrom1DArray(float oneDArray[], int x, int y, int z, int nx, int ny, int nz){
	int index = x*ny*nz + y*nz + z;
	return oneDArray[index];
}

float get2DValFrom1DArray(float oneDArray[], int x, int y, int nx, int ny){
	int index = x*ny + y;
	return oneDArray[index];
}

float tenRetrieveH(float* oneDArray, int nx, int ny, int nz, int x, int y, int z){
	int index = x*ny*nz + y*nz + z;
	return oneDArray[index];
}

void printSpectra(float *** integralMatrix, int nPixX, int nPixY){
    printf("\nX = [%f, ", integralMatrix[0][0][0]);
    for (int i=1; i<255; i++) {
        printf("%f, ", integralMatrix[0][0][i]);
    }
    printf("%f]\n ", integralMatrix[0][0][255]); 
    printf("\nY = [%f, ", integralMatrix[nPixX/6][nPixY/6][0]);
    for (int i=1; i<255; i++) {
        printf("%f, ", integralMatrix[nPixX/6][nPixY/6][i]);
    }
    printf("%f]\n ",integralMatrix[nPixX/6][nPixY/6][255]);
    printf("\nZ = [%f, ", integralMatrix[2*nPixX/6][2*nPixY/6][0]);
    for (int i=1; i<255; i++) {
        printf("%f, ", integralMatrix[2*nPixX/6][2*nPixY/6][i]);
    }
    printf("%f]\n ",integralMatrix[2*nPixX/6][2*nPixY/6][255]);
    printf("\nN = [%f, ", integralMatrix[3*nPixX/6][3*nPixY/6][0]);
    for (int i=1; i<255; i++) {
        printf("%f, ", integralMatrix[3*nPixX/6][3*nPixY/6][i]);
    }
    printf("%f]\n ",integralMatrix[3*nPixX/6][3*nPixY/6][255]);
    printf("\nM = [%f, ", integralMatrix[4*nPixX/6][4*nPixY/6][0]);
    for (int i=1; i<255; i++) {
        printf("%f, ", integralMatrix[4*nPixX/6][4*nPixY/6][i]);
    }
    printf("%f]\n ",integralMatrix[4*nPixX/6][4*nPixY/6][255]);
    printf("\nW = [%f, ", integralMatrix[5*nPixX/6][5*nPixY/6][0]);
    for (int i=1; i<255; i++) {
        printf("%f, ", integralMatrix[5*nPixX/6][5*nPixY/6][i]);
    }
    printf("%f]\n ",integralMatrix[5*nPixX/6][5*nPixY/6][255]);
    printf("\nT = [%f, ", integralMatrix[nPixX-1][nPixY-1][0]);
    for (int i=1; i<255; i++) {
        printf("%f, ", integralMatrix[nPixX-1][nPixY-1][i]);
    }
    printf("%f]\n ",integralMatrix[nPixX-1][nPixY-1][255]);
}

