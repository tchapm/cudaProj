//
//  cudaMethods.cpp
//  cudaProj
//
//  Created by Tyler Chapman on 2/21/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//
#include "cudaMethods.cuh"
#include <iostream>

__global__ void cudFunctionTest(float *testArr_d){
//    printf("test");
    int i=blockDim.x * blockIdx.x + threadIdx.x;
    testArr_d[i] = 1.0;
}

__device__ void modifyTensor(float* oneDArray, int index){
	oneDArray[index] = 100.0;
}

__device__ float tenRetrieveD(float* oneDArray, int nx, int ny, int nz, int x, int y, int z){
	int index = x*ny*nz + y*nz + z;
	return oneDArray[index];
}

__device__ float linearInterpZ(float* interpMat, int x, int y, float z, constants theConst){
	int z1, z2;
	int nx = theConst.nx;
	int ny = theConst.ny;
	int nz = theConst.nz;
	float ten1, ten2, interpValue;
	z1 = (int)z;
	z2 = z1+1;
	ten1 = __powf(10,tenRetrieveD(interpMat, nx, ny, nz, x, y, z1));
	ten2 = __powf(10,tenRetrieveD(interpMat, nx, ny, nz, x, y, z2));
	interpValue = ten1 + (ten2 -ten1)*(z - z1);
	return interpValue;
	
}

__device__ int getLowerIndex(float* axisArr, float value, int arrLength){
	for(int i=0; i<arrLength; i++){
		if (value<axisArr[i]) {
			return i-1;
		}
	}
	return arrLength-2;
}

__device__ float bilinInterpVal(float* interpMat, float y, float z, int energyBin, float* tempGrid, float* metalGrid, constants theConst){
	int y1, y2, z1, z2;
	float temp1, temp2, R1, R2, interpValue; 
	float ten11, ten12, ten21, ten22;
	float y1val, y2val, z1val, z2val;
	int nx = theConst.binCenterSize; //designed to interpolate the rebinned cooling function
	int ny = theConst.tGridSize;
	int nz = theConst.mGridSize;
	
	y1 = getLowerIndex(tempGrid, y, theConst.tGridSize);
	y2 = y1+1;
	z1 = getLowerIndex(metalGrid, z, theConst.mGridSize);
	z2 = z1+1;
	
	y1val = tempGrid[y1];
	y2val = tempGrid[y2];
	z1val = metalGrid[z1];
	z2val = metalGrid[z2];
	
	temp1 = (y2val-y)/(y2val-y1val);
	temp2 = (y-y1val)/(y2val-y1val);
	ten11 = __powf(10,tenRetrieveD(interpMat, nx, ny, nz, energyBin, y1, z1));//tenRetrieveD(interpMat, nx, ny, nz, energyBin, y1, z1);//
	ten21 = __powf(10,tenRetrieveD(interpMat, nx, ny, nz, energyBin, y2, z1));//tenRetrieveD(interpMat, nx, ny, nz, energyBin, y2, z1);//
	ten12 = __powf(10,tenRetrieveD(interpMat, nx, ny, nz, energyBin, y1, z2));//tenRetrieveD(interpMat, nx, ny, nz, energyBin, y1, z2);//
	ten22 = __powf(10,tenRetrieveD(interpMat, nx, ny, nz, energyBin, y2, z2));//tenRetrieveD(interpMat, nx, ny, nz, energyBin, y2, z2);
	R1 = temp1*ten11 + temp2*ten21;
	R2 = temp1*ten12 + temp2*ten22;
	interpValue = ((z2val-z)/(z2val-z1val))*R1 + ((z-z1val)/(z2val-z1val))*R2;
	
	return interpValue;
}


__global__ void integrate(float* rebinCool, float* tempArr, float* metalArr, float* emmArr, float* integral, float* tempGrid, float* metalGrid, constants theConst, bool debugging){
	//integrate from depth/2 to -depth/2
	
//    printf("!!!!!!!!!!!!!!!!!!!test!!!!!!!!!!!!!!!!!");
    int i=blockDim.x * blockIdx.x + threadIdx.x; //threadIdx.x is the channel/energy-bin
    int j= blockIdx.x;   
	int x, y;
	x = j/theConst.nPixX;
	y = j%theConst.nPixY;
	int energyBin = threadIdx.x;
	int a, b, n, step=1;
	float tFunct, h, actErr, last, nextVal;
	float T, Z, rebinA, rebinB;
	float prevStep[200];
	//b = theConst.depth/2 + theConst.nz/2;
	//a = theConst.nz/2 - theConst.depth/2;
	b = theConst.depth + theConst.nz/2-1;
	a = theConst.nz/2-1;
	h = b-a;
	actErr = 1.0;
	last = 1.E30;
	nextVal = 0.0;
	n = 1;
	T = __powf(10,tenRetrieveD(tempArr, theConst.nx, theConst.ny, theConst.nz, x, y, b));
	Z = __powf(10,tenRetrieveD(metalArr, theConst.nx, theConst.ny, theConst.nz, x, y, b));	
	rebinB = bilinInterpVal(rebinCool, T, Z, energyBin, tempGrid, metalGrid, theConst);
	T = __powf(10,tenRetrieveD(tempArr, theConst.nx, theConst.ny, theConst.nz, x, y, a));
	Z = __powf(10,tenRetrieveD(metalArr, theConst.nx, theConst.ny, theConst.nz, x, y, a));
	rebinA = bilinInterpVal(rebinCool, T, Z, energyBin, tempGrid, metalGrid, theConst);
	tFunct = 0.5*__powf(__powf(10,tenRetrieveD(emmArr, theConst.nx, theConst.ny, theConst.nz, x, y, a)), 2.0)* rebinA;
	tFunct += 0.5*__powf(__powf(10,tenRetrieveD(emmArr, theConst.nx, theConst.ny, theConst.nz, x, y, b)), 2.0)*rebinB;
    //iterate until convergence of integral
    prevStep[0] = tFunct;
	while (actErr>=0.1) {
		step=step*2;
		h = float(b-a)/step;
		nextVal = 0.0;
		for (int l=1; l<step; l=l+2) {
			T = linearInterpZ(tempArr, x, y, l*h+a, theConst);
            Z = linearInterpZ(metalArr, x, y, l*h+a, theConst);
			nextVal+=bilinInterpVal(rebinCool, T, Z, energyBin, tempGrid, metalGrid, theConst)*__powf(__powf(10,tenRetrieveD(emmArr, theConst.nx, theConst.ny, theConst.nz, x, y, l*h+a)), 2.0);
		}
		nextVal+=prevStep[n-1];
		prevStep[n]=nextVal;
		nextVal=h*(nextVal);
		actErr=fabs(last-nextVal);
		last=nextVal;
		n++;
	}
    
    if (debugging) {
        integral[0]=T;
        integral[1]=energyBin;
        integral[2]=__powf(__powf(10,tenRetrieveD(emmArr, theConst.nx, theConst.ny, theConst.nz, x, y, 2*h+a)), 2.0);
        integral[3]=Z;
        integral[4]=bilinInterpVal(rebinCool, T, Z, energyBin, tempGrid, metalGrid, theConst);
        integral[5]=tFunct;
        integral[6]=__powf(__powf(10,tenRetrieveD(emmArr, theConst.nx, theConst.ny, theConst.nz, x, y, b)), 2.0);
        integral[7] = rebinA;
        integral[8] = n;
        integral[9] = last;
        integral[10] = prevStep[1];
    }else {
    ////place integrations into the 1D array by thread number
        integral[i] = last;
    }	 
}

__global__ void tempInit(void(*f)(float), float * temperature, constants theConst){
	/*int i=blockDim.x * blockIdx.x + threadIdx.x; //threadIdx.x is the channel or energy-bin	
	int j= blockIdx.x;   
	int x, y;
	x = j/theConst.nx;
	y = j%theConst.ny;
    int z = threadIdx.x;
    float l, A, B, C, rNot;
    A = 1.0;
    B = 1.0;
    C = 1.0;
    float a, b, actErr, h;
    int n, step;
    //    a = 1.0;
    b = 10.0; //need to check conversion .01 pix/Mpc
    //    h=b-a;
    rNot = 0.5; //500 kPc * 0.01 pixel/Mpc
    float tFunct, nextVal, last;
    float prevStep[200];
	
    n=1;
    nextVal=0.0;
    last = 1.E30;
    actErr = 1.0;
    step = 1;
    l = 1.0;//0.01*__powf(__powf(float(x-(theConst.nx/2-1))/A,2)+__powf(float(y-(theConst.ny/2-1))/B,2)+__powf(float(z-(theConst.nz/2-1))/C,2),0.5);
    tFunct = l;//0.5*(*f)(rNot,l)
    (*f)(rNot);
    //    tFunct += 0.5*(*f)(rNot,b);
    //    prevStep[0] = tFunct;
    //    while (actErr>=0.01) {
    //        step=step*2;
    //        nextVal = 0.0;
    //        h = float(b-l)/step;		
    //        for (int k=1; k<step; k=k+2) {
    //            nextVal += (*f)(rNot,k*h+l);
    //        }
    //        nextVal+=prevStep[n-1];
    //        prevStep[n]=nextVal;
    //        nextVal=h*(nextVal);
    //        actErr=fabs(last-nextVal);
    //        last=nextVal;
    //        n++;
    //    }
    temperature[i]= z;
    */
}

__device__ float function (float rNot, float l) {
    return 1.0;//__powf(rNot+l,-5.0);
}

__device__ void function2 (float x) {
    1.0;
}

