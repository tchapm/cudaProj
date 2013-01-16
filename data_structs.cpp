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
#include "Structs.h"
//#include <sciutils.h>
//#include <jaco.h>

using namespace std;

double **sci_fmatrixD(int nx, int ny)
{
	double **x=NULL;
	long i;
	
	if (nx > 0 && ny > 0) 
		x = (double **)malloc(nx*sizeof(double *));
	if (x == NULL){ 
		printf("Could not allocate float matrix with %d rows",nx);
		perror("errno");
		exit(1);
	}
	
	x[0] = (double *)malloc(nx*ny*sizeof(double));
	if (x[0] == NULL){ 
		printf("Could not allocate %d%d float matrix",nx,ny);
		perror("errno");
		exit(1);
		
	}
	for (i = 1; i < nx; ++i) 
		x[i] = x[0]+i*ny;
	
	return x;
}
//needed to be float for cooling function calculation
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
		printf("Could not allocate double matrix with %d rows",nx);
		perror("errno");
		exit(1);
	}
	
	x[0][0] = (double *)malloc(nx*ny*nz*sizeof(double));
	if (x[0] == NULL){ 
		printf("Could not allocate %d%d double matrix",nx,ny);
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

double *tensorTo1DArray(double ***tensor3D, int nx, int ny, int nz){
	double *oneDArray = new double[nx*ny*nz];
	
	int arrIndex = 0;
	for (int i=0; i<nx; i++) {
		for (int j=0; j<ny; j++) {
			for (int k=0; k<nz; k++) {
				oneDArray[arrIndex]=(double)tensor3D[i][j][k];
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
			arrIndex++;
		}
	}
	return oneDArray;
}


double *metalArrInit(int nEll, double a_ell, double b_ell, double ellMax){
	double *metalArr = new double[nEll];
    double *ellArr = ellArrInit(nEll, ellMax);
    double offset;
    for(int k=0; k<nEll; k++){
        //why 0.907????
        offset = k*0.907;
        metalArr[k] = log10((nEll-offset)/nEll);
//        printf("metalArr[%d] = %f\n", k,  powf(10,metalArr[k]));
    }
	return metalArr;
}

double *emmArrInit(int nEll, double a_ell, double b_ell, double ellMax){
	double *emmArr = new double[nEll];
    double rNot = 0.5; //500 kPc * 0.01 pixel/Mpc
    double *ellArr = new double[nEll];
    ellArr = ellArrInit(nEll, ellMax);
    for(int k=0; k<nEll; k++){
        emmArr[k] = log10(pow(rNot+ellArr[k],-3));
        //        printf("emmArr[%d] = %f\n", k,  emmArr[k]);
    }
	return emmArr;
}

double* tempArrInit(int nEll, double a_ell, double b_ell, double ellMax){
	double *tempArr = new double[nEll];
    double rNot = 0.5; //500 kPc * 0.01 pixel/Mpc
    double *ellArr = new double[nEll];
    ellArr = ellArrInit(nEll, ellMax);
    for(int k=0; k<nEll; k++){
        tempArr[k]=log10(6.7*pow(rNot+ellArr[k],-1.0));
        //        printf("tempArr[%d] = %f\n", k,  tempArr[k]);
    }
    
	return tempArr;
}

//initializers not using ellArr
/*
double *emmArrInit(int nx, int ny, int nz){
	int x, y, z;
	double r;
	double *emmArr = new double[nx*ny*nz];
    double l, A, B, C, rNot;
    A = 1.0;
    B = 1.0;
    C = 1.0;
    rNot = 0.5;
	for (x=0; x<nx; x++) {
		for(y=0; y<ny; y++){
			for (z=0; z<nz; z++) {
                l =powf((powf(double(x-(nx/2-1))/A,2)+powf(double(y-(ny/2-1))/B,2)+powf(double(z-(nz/2-1))/C,2)), 0.5);
				r=powf((powf(x-(nx/2-1),2)+powf(y-(ny/2-1),2)+powf(z-(nz/2-1), 2)), 0.5);
				emmArr[x*ny*nz + y*nz + z] = log10(powf(rNot+0.01*l,-3));
			}
		}
	}
	return emmArr;
}

double *metalArrInit(int nx, int ny, int nz){
	int x, y, z;
	double r;
	double *metalArr = new double[nx*ny*nz];
    
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


double *tempArrInit(int nx, int ny, int nz){
	int x, y, z;
	double r;
	double *tempArr = new double[nx*ny*nz];
    double l, A, B, C, rNot;
    A = 1.0;
    B = 1.0;
    C = 1.0;
    rNot = 0.5; //500 kPc * 0.01 pixel/Mpc
	for (x=0; x<nx; x++) {
		for(y=0; y<ny; y++){
			for (z=0; z<nz; z++) {
                l =powf((powf(double(x-(nx/2-1))/A,2)+powf(double(y-(ny/2-1))/B,2)+powf(double(z-(nz/2-1))/C,2)), 0.5);
				r=powf((powf(x-(nx/2-1),2)+powf(y-(ny/2-1),2)+powf(z-(nz/2-1), 2)), 0.5);
				tempArr[x*ny*nz + y*nz + z]=log10(6.7*powf(rNot+0.01*l,-1.0));
			}
		}
	}
	return tempArr;
}

*/

double* energyArrInit(int rebinSize){
	int x;
	double eMin = 0.2; //in Kev
	double eMax = 10.0;
	double binWidth = (eMax-eMin)/rebinSize;
	double *energyArr = new double[rebinSize];
	energyArr[0] = eMin;
	for (x=1; x<rebinSize; x++) {
		energyArr[x] = energyArr[x-1]+binWidth;
	}
	return energyArr;
	
}

double* ellArrInit(int nEll, double ellMax){
	int x;
	double eMin = 0.1; 
	double eMax = ellMax;
	double binWidth = (eMax-eMin)/nEll;
	double *ellArr = new double[nEll];
	ellArr[0] = binWidth;
	for (x=1; x<nEll; x++) {
		ellArr[x] = ellArr[x-1]+binWidth;
	}	
//    printf("ell range: %f to %f over %d bins\n", eMin, ellArr[nEll-1], nEll);
	return ellArr;
}


double* tempGridInit(int tGridSize, double *tempAxis){
    double *tempGrid_h = (double*) calloc(tGridSize, sizeof(double));
    for(int i =0; i<tGridSize; i++){
		tempGrid_h[i] = pow(10,tempAxis[i]);
	}
    return tempGrid_h;
}

double* metalGridInit(int mGridSize, double *metalAxis){
    double *metalGrid_h = (double*) calloc(mGridSize, sizeof(double));
    for(int i =0; i<mGridSize; i++){
		metalGrid_h[i] = pow(10,metalAxis[i]);
	}
    return metalGrid_h;
}

double* rotMatInit(double theta, double phi, double psi){
    double *rotMat = new double[9];
    rotMat[0] = cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi);
    rotMat[1] = cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi);
    rotMat[2] = sin(psi)*sin(theta);
    rotMat[3] = -sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi);
    rotMat[4] = -sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi);
    rotMat[5] = cos(psi)*sin(theta);
    rotMat[6] = sin(theta)*sin(phi);
    rotMat[7] = -sin(theta)*cos(phi);
    rotMat[8] = cos(theta);
    return rotMat;
}

double*** makeIntegralMatrix(double *integral_h, int nPixX, int nPixY, int binCenterSize){
    double ***integralMatrix;
    integralMatrix = sci_ftensorD(nPixX, nPixY, binCenterSize);
    for (int x=0; x<nPixX; x++) {
        for (int y=0; y<nPixY; y++) {
            for (int z=0; z<binCenterSize; z++) {
                integralMatrix[x][y][z] = log10(tenRetrieveH(integral_h, nPixX, nPixY, binCenterSize, x, y, z));
            }
        }
    }
    return integralMatrix;
}

double get3DValFrom1DArray(double oneDArray[], int x, int y, int z, int nx, int ny, int nz){
	int index = x*ny*nz + y*nz + z;
	return oneDArray[index];
}

double get2DValFrom1DArray(double oneDArray[], int x, int y, int nx, int ny){
	int index = x*ny + y;
	return oneDArray[index];
}

double tenRetrieveH(double* oneDArray, int nx, int ny, int nz, int x, int y, int z){
	int index = x*ny*nz + y*nz + z;
	return oneDArray[index];
}

void printSpectra(double*** integralMatrix, double* energyArr, constants theConst){
    printf("\nX = [%f, ", integralMatrix[0][0][0]);
    for (int i=1; i<theConst.binCenterSize-1; i++) {
        printf("%f, ", integralMatrix[0][0][i]);
    }
    printf("%f]\n ", integralMatrix[0][0][theConst.binCenterSize-1]); 
    printf("\nY = [%f, ", integralMatrix[theConst.nPixX/6][theConst.nPixY/6][0]);
    for (int i=1; i<theConst.binCenterSize-1; i++) {
        printf("%f, ", integralMatrix[theConst.nPixX/6][theConst.nPixY/6][i]);
    }
    printf("%f]\n ",integralMatrix[theConst.nPixX/6][theConst.nPixY/6][theConst.binCenterSize-1]);
    printf("\nZ = [%f, ", integralMatrix[2*theConst.nPixX/6][2*theConst.nPixY/6][0]);
    for (int i=1; i<theConst.binCenterSize-1; i++) {
        printf("%f, ", integralMatrix[2*theConst.nPixX/6][2*theConst.nPixY/6][i]);
    }
    printf("%f]\n ",integralMatrix[2*theConst.nPixX/6][2*theConst.nPixY/6][theConst.binCenterSize-1]);
    printf("\nN = [%f, ", integralMatrix[3*theConst.nPixX/6][3*theConst.nPixY/6][0]);
    for (int i=1; i<theConst.binCenterSize-1; i++) {
        printf("%f, ", integralMatrix[3*theConst.nPixX/6][3*theConst.nPixY/6][i]);
    }
    printf("%f]\n ",integralMatrix[3*theConst.nPixX/6][3*theConst.nPixY/6][theConst.binCenterSize-1]);
    printf("\nM = [%f, ", integralMatrix[4*theConst.nPixX/6][4*theConst.nPixY/6][0]);
    for (int i=1; i<theConst.binCenterSize-1; i++) {
        printf("%f, ", integralMatrix[4*theConst.nPixX/6][4*theConst.nPixY/6][i]);
    }
    printf("%f]\n ",integralMatrix[4*theConst.nPixX/6][4*theConst.nPixY/6][theConst.binCenterSize-1]);
    printf("\nW = [%f, ", integralMatrix[5*theConst.nPixX/6][5*theConst.nPixY/6][0]);
    for (int i=1; i<theConst.binCenterSize-1; i++) {
        printf("%f, ", integralMatrix[5*theConst.nPixX/6][5*theConst.nPixY/6][i]);
    }
    printf("%f]\n ",integralMatrix[5*theConst.nPixX/6][5*theConst.nPixY/6][theConst.binCenterSize-1]);
    printf("\nT = [%f, ", integralMatrix[theConst.nPixX-1][theConst.nPixY-1][0]);
    for (int i=1; i<theConst.binCenterSize-1; i++) {
        printf("%f, ", integralMatrix[theConst.nPixX-1][theConst.nPixY-1][i]);
    }
    printf("%f]\n ",integralMatrix[theConst.nPixX-1][theConst.nPixY-1][theConst.binCenterSize-1]);
    printf("\nenergyBin = [%f, ", energyArr[0]);
    for (int i=1; i<theConst.binCenterSize-1; i++) {
        printf("%f, ", energyArr[i]);
    }
    printf("%f]\n ",energyArr[theConst.binCenterSize-1]);
    double pixToMpc = 0.01;
    printf("XMpc = %f\n", pixToMpc*0*theConst.nPixY/6*pow(2, 0.5));
    printf("YMpc = %f\n", pixToMpc*1*theConst.nPixY/6*pow(2, 0.5));
    printf("ZMpc = %f\n", pixToMpc*2*theConst.nPixY/6*pow(2, 0.5));
    printf("NMpc = %f\n", pixToMpc*3*theConst.nPixY/6*pow(2, 0.5));
    printf("MMpc = %f\n", pixToMpc*4*theConst.nPixY/6*pow(2, 0.5));
    printf("WMpc = %f\n", pixToMpc*5*theConst.nPixY/6*pow(2, 0.5));
    printf("TMpc = %f\n", pixToMpc*6*theConst.nPixY/6*pow(2, 0.5));
}


void sumSpectra(double*** integralMatrix, double* energyArr, constants theConst){
    double sum = 0.0;
    double minEnergy = 0.2, maxEnergy = 10.0;
    double eBin = (maxEnergy-minEnergy)/theConst.binCenterSize;
    double pixToMpc = 0.01;
    for (int i=0; i<theConst.binCenterSize; i++) {
        sum += pow(10,integralMatrix[0*theConst.nPixX/128][0*theConst.nPixY/128][i]);
    }
    printf("\nFlux = [%f ", sum*eBin);
    for (int j=1; j<128; j++) {
        sum = 0.0;
        for (int i=0; i<theConst.binCenterSize; i++) {
            sum += pow(10,integralMatrix[j*theConst.nPixX/128][j*theConst.nPixY/128][i]);
        }
        printf(", %f", sum*eBin);
    }
    printf("]");
    printf("\nposNum = [%d ", 0);
    for (int i=1; i<128; i++) {
        printf(", %f", i*pixToMpc*pow(2,0.5));
    }
    printf("]\n");
    
}

void plotImage(double*** integralMatrix, double* energyArr, constants theConst, char* fileName){
    double sum = 0.0;
    double minEnergy = 0.2, maxEnergy = 10.0;
    double eBin = (maxEnergy-minEnergy)/theConst.binCenterSize;
    FILE *file; 
    //open file to append
    file = fopen(fileName ,"w"); 
    fprintf(file,"%d\n", theConst.nPixX*theConst.nPixY);
//    printf("\nFlux = [%f ", sum*eBin);
    int i=0;
    for (int k=0; k<theConst.nPixX; k++) {
        for (int j=0; j<theConst.nPixY; j++) {
            sum = 0.0;
            for (i=0; i<theConst.binCenterSize; i++) {
                sum += pow(10,integralMatrix[k][j][i]);
            }
            fprintf(file," %f", sum*eBin);
//            printf(", %f", sum*eBin);
        }
    }
//    fprintf(file,"%s","]\n");
//    printf("]\n");
    fclose(file); 
}
