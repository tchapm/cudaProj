//
//  cudaMethods.cpp
//  cudaProj
//
//  Created by Tyler Chapman on 2/21/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//
#include "cudaMethods.cuh"
#include <iostream>

__global__ void cudFunctionTest(double *testArr_d){
    int i=blockDim.x * blockIdx.x + threadIdx.x;
    testArr_d[i] = 1.0;
}

__device__ double tenRetrieve(double* oneDArray, int nx, int ny, int nz, int x, int y, int z){
	int index = x*ny*nz + y*nz + z;
	return oneDArray[index];
}

//returns log10
__device__ double linearInterpEll(double* interpMat, double ell, constants theConst){
	int ell1, ell2;
	double ten1, ten2, interpValue;
    //find raw value on index array
    double ellIndex = ell*theConst.n_ell/theConst.ellMax;
	ell1 = (int)(ellIndex);
	ell2 = ell1+1;
    //retrieve the array value for two indexes
	ten1 = __powf(10,interpMat[ell1]);
	ten2 = __powf(10,interpMat[ell2]);
	interpValue = ten1 + (ten1 - ten2)*(ellIndex - (int)ellIndex);
	return log10(interpValue);	
}

__device__ int getLowerIndex(double* axisArr, double value, int arrLength){
	for(int i=0; i<arrLength; i++){
		if (value<axisArr[i]) {
			return i-1;
		}
	}
	return arrLength-2;
}

//returns unLog10
__device__ double bilinInterpVal(double* interpMat, double y, double z, int energyBin, double* tempGrid, double* metalGrid, constants theConst){
	int y1, y2, z1, z2;
	double temp1, temp2, R1, R2, interpValue; 
	double ten11, ten12, ten21, ten22;
	double y1val, y2val, z1val, z2val;
	int nx = theConst.binCenterSize; 
	int ny = theConst.tGridSize;
	int nz = theConst.mGridSize;
	
	y1 = getLowerIndex(tempGrid, __powf(10,y), theConst.tGridSize);
	y2 = y1+1;
	z1 = getLowerIndex(metalGrid, __powf(10,z), theConst.mGridSize);
	z2 = z1+1;

	y1val = tempGrid[y1];
	y2val = tempGrid[y2];
	z1val = metalGrid[z1];
	z2val = metalGrid[z2];
	y = __powf(10,y);
    z = __powf(10,z);
	temp1 = (y2val-y)/(y2val-y1val);
	temp2 = (y-y1val)/(y2val-y1val);
	ten11 = __powf(10,tenRetrieve(interpMat, nx, ny, nz, energyBin, y1, z1));
	ten21 = __powf(10,tenRetrieve(interpMat, nx, ny, nz, energyBin, y2, z1));
	ten12 = __powf(10,tenRetrieve(interpMat, nx, ny, nz, energyBin, y1, z2));
	ten22 = __powf(10,tenRetrieve(interpMat, nx, ny, nz, energyBin, y2, z2));
	R1 = temp1*ten11 + temp2*ten21;
	R2 = temp1*ten12 + temp2*ten22;
	interpValue = ((z2val-z)/(z2val-z1val))*R1 + ((z-z1val)/(z2val-z1val))*R2;
	
	return interpValue;
}

__device__ double getEll(double x, double y, double z, double* rotMat, double a, double b){
    double xPrime, yPrime, zPrime;
    xPrime = rotMat[0]*x + rotMat[1]*y + rotMat[2]*z;
    yPrime = rotMat[3]*x + rotMat[4]*y + rotMat[5]*z;
    zPrime = rotMat[6]*x + rotMat[7]*y + rotMat[8]*z;
    return __powf(xPrime*xPrime + yPrime*yPrime/(a*a) + zPrime*zPrime/(b*b),0.5);
}


__global__ void integrate(double* rebinCool, double* tempArr, double* metalArr, double* emmArr, double* integral, double* tempGrid, double* metalGrid, double* rotMat, constants theConst, bool debugging){
	
    int i=blockDim.x * blockIdx.x + threadIdx.x; //threadIdx.x is the channel/energy-bin
    int j= blockIdx.x;   
	int x, y;
    //integrate from depth/2 to -depth/2
    x = j/theConst.nPixX-(theConst.nPixX/2);
	y = j%theConst.nPixY-(theConst.nPixY/2);
    double pixToMpc = 0.01;

    double xMpc = pixToMpc*x;
    double yMpc = pixToMpc*y;
	int energyBin = threadIdx.x;
	int n;
    double a, b;
    double step=1.0;
	double tFunct, h, actErr, last, nextVal;
	double T, Z, rebinA, rebinB;
	double prevStep[200];
    double ell;
	b = theConst.nz/2;
	a = -theConst.nz/2;
	h = b-a;
	actErr = 1.0;
	last = 1.E30;
	nextVal = 0.0;
	n = 1;
    
    ell = getEll(xMpc, yMpc, a, rotMat, theConst.a_ell, theConst.b_ell);
    if (ell>theConst.ellMax) {
        integral[i] = 1.0;
    }
    T = linearInterpEll(tempArr, ell, theConst);
	Z = linearInterpEll(metalArr, ell, theConst);
    double ell1 = getEll(xMpc, yMpc, a, rotMat, theConst.a_ell, theConst.b_ell);
    double T1 = linearInterpEll(tempArr, ell, theConst);
    double Z1 = linearInterpEll(metalArr, ell, theConst);
	rebinA = bilinInterpVal(rebinCool, T, Z, energyBin, tempGrid, metalGrid, theConst);
	tFunct = 0.5*__powf(__powf(10,linearInterpEll(emmArr, ell, theConst)), 2.0)*rebinA;
    ell = getEll(xMpc, yMpc, b, rotMat, theConst.a_ell, theConst.b_ell);
    if (ell>theConst.ellMax) {
        integral[i] = 1.0;
    }
	T = linearInterpEll(tempArr, ell, theConst);
	Z = linearInterpEll(metalArr, ell, theConst);
	rebinB = bilinInterpVal(rebinCool, T, Z, energyBin, tempGrid, metalGrid, theConst);
	tFunct += 0.5*__powf(__powf(10,linearInterpEll(emmArr, ell, theConst)), 2.0)*rebinB;
    //iterate until convergence of integral
    prevStep[0] = tFunct;
	while (actErr>=0.1) {
		step=step*2.0;
		h = double(b-a)/step;
		nextVal = 0.0;
		for (int l=1; l<step; l=l+2) {
            ell = getEll(xMpc, yMpc, l*h+a, rotMat, theConst.a_ell, theConst.b_ell);
            if (ell>theConst.ellMax) {
                integral[i] = 1.0;
                break;
            }
			T = linearInterpEll(tempArr, ell, theConst);
            Z = linearInterpEll(metalArr, ell, theConst);
			nextVal+=bilinInterpVal(rebinCool, T, Z, energyBin, tempGrid, metalGrid, theConst)*__powf(__powf(10,linearInterpEll(emmArr, ell, theConst)), 2.0);
		}
		nextVal+=prevStep[n-1];
		prevStep[n]=nextVal;
		nextVal=h*(nextVal);
		actErr=fabs(last-nextVal);
		last=nextVal;
		n++;
        if (debugging && n>4) {
            break;
        }
	}
    
    if (debugging) {
        integral[0]=xMpc;
        integral[1]=energyBin;
        integral[2]=ell1;
        integral[3]=Z1;
        integral[4]=linearInterpEll(tempArr, ell, theConst);
        integral[5]=tFunct;
        integral[6]=linearInterpEll(emmArr, ell, theConst);
        integral[7] = rebinA;
        integral[8] = n;
        integral[9] = last;
        integral[10] = j;
        integral[11] = a;
        integral[12] = b;
        integral[13] = T1;
        

    }else {
        //place integrations into the 1D array by thread number
        //if outside the grid then make minescule 
        if (ell>theConst.ellMax || integral[i] == 1.0) {
            integral[i] = 1.0E-10;
        }else{
            integral[i] = last;
        }
    }	 
}


