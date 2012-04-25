
/*
 *  cud3Dsim.cpp
 *  3Dsim
 *
 *  Created by Tyler Chapman on 10/16/10.
 *
 */
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <cutil_inline.h>
#include <iomanip>
#include <errno.h>
#include <vector>
#include <stdlib.h>
#include <fcntl.h>
#include "cooling_functions.cpp"
#include "Structs.h"
#include "cuda_data_structs.cu"
#include "cudaMethods.cu"
#include "jaco.h"

using namespace std;

double *metalArr_h;
double *metalArr_d;
double *emmArr_h;
double *emmArr_d;
double *tempArr_h;
double *tempArr_d;
double *energyArr_h;
//double *energyArr_d;
double *rebinArr_h;
double *rebinArr_d;
double *integral_h;
double *integral_d;
double *metalGrid_h;
double *metalGrid_d;
double *tempGrid_h;
double *tempGrid_d;
double *rotMat_h;
double *rotMat_d;


void Cleanup(void)
{
    // Free device memory
	if (metalArr_d) {
		cudaFree(metalArr_d);
	}if (emmArr_d) {
		cudaFree(emmArr_d);
	}if (tempArr_d) {
		cudaFree(tempArr_d);
	}if (rebinArr_d) {
		cudaFree(rebinArr_d);
	}if (integral_d) {
		cudaFree(integral_d);
	}if (metalGrid_d) {
		cudaFree(metalGrid_d);
	}if (tempGrid_d) {
		cudaFree(tempGrid_d);
	}if (rotMat_d) {
		cudaFree(rotMat_d);
	}
	
	// Free host memory
	if (metalArr_h) {
		free(metalArr_h);
	}if (emmArr_h) {
		free(emmArr_h);
	}if (tempArr_h) {
		free(tempArr_h);
	}if (energyArr_h) {
		free(energyArr_h);
	}if (rebinArr_h) {
		free(rebinArr_h);
	}if (integral_h) {
		free(integral_h);
	}if (metalGrid_h) {
		free(metalGrid_h);
	}if (tempGrid_h) {
		free(tempGrid_h);
	}if (rotMat_h) {
		free(rotMat_h);
	}
	
    cutilSafeCall( cudaThreadExit() );	
    exit(0);
}


//
//__global__ void integrate(float* rebinCool, float* tempArr, float* metalArr, float* emmArr, float* integral, float* tempGrid, float* metalGrid, constants theConst, bool debugging){
//	//integrate from depth/2 to -depth/2
//	int i=blockDim.x * blockIdx.x + threadIdx.x; //threadIdx.x is the channel/energy-bin	
//	int j= blockIdx.x;   
//	int x, y;
//	x = j/theConst.nPixX;
//	y = j%theConst.nPixY;
//	int energyBin = threadIdx.x;
//	int a, b, n, step=1;
//	float tFunct, h, actErr, last, nextVal;
//	float T, Z, rebinA, rebinB;
//	float prevStep[200];
//	//b = theConst.depth/2 + theConst.nz/2;
//	//a = theConst.nz/2 - theConst.depth/2;
//	b = theConst.depth + theConst.nz/2-1;
//	a = theConst.nz/2-1;
//	h = b-a;
//	actErr = 1.0;
//	last = 1.E30;
//	nextVal = 0.0;
//	n = 1;
//    
//	T = __powf(10,tenRetrieveD(tempArr, theConst.nx, theConst.ny, theConst.nz, x, y, b));
//	Z = __powf(10,tenRetrieveD(metalArr, theConst.nx, theConst.ny, theConst.nz, x, y, b));	
//	rebinB = bilinInterpVal(rebinCool, T, Z, energyBin, tempGrid, metalGrid, theConst);
//	T = __powf(10,tenRetrieveD(tempArr, theConst.nx, theConst.ny, theConst.nz, x, y, a));
//	Z = __powf(10,tenRetrieveD(metalArr, theConst.nx, theConst.ny, theConst.nz, x, y, a));
//	rebinA = bilinInterpVal(rebinCool, T, Z, energyBin, tempGrid, metalGrid, theConst);
//	tFunct = 0.5*__powf(__powf(10,tenRetrieveD(emmArr, theConst.nx, theConst.ny, theConst.nz, x, y, a)), 2.0)* rebinA;
//	tFunct += 0.5*__powf(__powf(10,tenRetrieveD(emmArr, theConst.nx, theConst.ny, theConst.nz, x, y, b)), 2.0)*rebinB;
//    //iterate until convergence of integral
//    prevStep[0] = tFunct;
//	while (actErr>=0.1) {
//		step=step*2;
//		h = float(b-a)/step;
//		nextVal = 0.0;
//		for (int l=1; l<step; l=l+2) {
//			T = linearInterpZ(tempArr, x, y, l*h+a, theConst);
//            Z = linearInterpZ(metalArr, x, y, l*h+a, theConst);
//			nextVal+=bilinInterpVal(rebinCool, T, Z, energyBin, tempGrid, metalGrid, theConst)*__powf(__powf(10,tenRetrieveD(emmArr, theConst.nx, theConst.ny, theConst.nz, x, y, l*h+a)), 2.0);
//		}
//		nextVal+=prevStep[n-1];
//		prevStep[n]=nextVal;
//		nextVal=h*(nextVal);
//		actErr=fabs(last-nextVal);
//		last=nextVal;
//		n++;
//	}
//    
//    if (debugging) {
//        integral[0]=T;
//        integral[1]=energyBin;
//        integral[2]=i;
//        integral[3]=Z;
//        integral[4]=j;
//        integral[5]=tFunct;
//        integral[6]=__powf(__powf(10,tenRetrieveD(emmArr, theConst.nx, theConst.ny, theConst.nz, x, y, b)), 2.0);
//        integral[7] = rebinA;
//        integral[8] = n;
//        integral[9] = last;
//    }
//    ////place integrations into the 1D array by thread number
//	integral[i] = last;	 
//}



double* sumArea(int x1, int x2, int y1, int y2, double*** completeArr, constants theConst){
    double* combinedRegion;
    combinedRegion = (double*) calloc(theConst.binCenterSize, sizeof(double));
    for (int k=0; k<theConst.binCenterSize; k++) {
        for (int i=x1; i<=x2; i++) {
            for (int j=y1; j<=y2; j++) {
                combinedRegion[k]+= powf(10,completeArr[i][j][k]);
            }
        }
        combinedRegion[k]=log10(combinedRegion[k]);
    }
    return combinedRegion;
}
void constInit(constants &theConst, jaco_state ja){
    theConst.eGridSize = ja.egridsize;
	theConst.tGridSize = ja.tgridsize;
	theConst.mGridSize = ja.mgridsize;
    theConst.tempAxis = ja.tempaxis;
    theConst.metalAxis = ja.metalaxis;
	theConst.depth = 2.0 * ja.rshock; //in Mpc
	theConst.nPixX = ja.nx;//128;
    theConst.nPixY = ja.ny;//128;
	theConst.pixToMpc = 1/(ja.angdist*ja.pixscale);//0.01;
	theConst.binCenterSize = ja.nlastbin;//256;
    theConst.theta = ja.theta;
    theConst.phi = ja.phi;
    theConst.epsilon = ja.epsilon;
    theConst.a_ell = ja.a_ell;
    theConst.b_ell = ja.b_ell;
    theConst.n_ell = ja.n_ell;
    theConst.rMax = ja.rshock;
    theConst.rebinnedCooling = ja.rebinnedcooling;
    theConst.metalArr=ja.Zprof;
    theConst.tempArr=ja.Tprof;
    theConst.emmArr=ja.nenhprof;
	theConst.accuracy=0.0001;
    //	theConst.nGrid=1;
//	theConst.nx=256;
//	theConst.ny=256;
//	theConst.nz=256;	
    
}


double *runSimulation(jaco_state ja, int x1, int x2, int y1, int y2){
    ////////////set constants and sizes of data structures
    clock_t sTime, eTime, sLoadTime, eLoadTime, sGPUTime, eGPUTime, sCoolTime, eCoolTime, sIntTime, eIntTime;
    /////////////start clock
    sTime = clock();
	constants theConst;
	constInit(theConst, ja);
    
//    file_reader coolingFile;
//	double centBin[theConst.binCenterSize];
//	strcpy(coolingFile.spectralCode, "/home/tchap/mekal.bin");
//	coolingFile.debug = 2;
//	read_cooling_function(coolingFile);
//    makeReBin(centBin, theConst.binCenterSize);
//	rebincoolingfunction(centBin, theConst.binCenterSize, coolingFile);
//    theConst.tGridSize = coolingFile.tGridSize;
//    theConst.mGridSize = coolingFile.mGridSize;
//    theConst.eGridSize = coolingFile.eGridSize;
//    theConst.tempAxis = coolingFile.tempAxis;
//    theConst.metalAxis = coolingFile.metalAxis;
    
	size_t sizeT = theConst.n_ell*sizeof(double);
	size_t sizeEn = theConst.binCenterSize*sizeof(double);
	size_t sizeInt = theConst.nPixX*theConst.nPixY*theConst.binCenterSize*sizeof(double);
    
	///////////cooling function variables

    size_t sizeTGrid = theConst.tGridSize*sizeof(double);
    size_t sizeMGrid = theConst.mGridSize*sizeof(double);
	
	/////////////create rebin and put into vectors to transfer to device	
    sCoolTime = clock();

    eCoolTime = clock();
    
    /////////////make host data structures
	printf("Loading input data tensors\n");
    sLoadTime = clock();
	rebinArr_h = tensorTo1DArray(theConst.rebinnedCooling, theConst.binCenterSize, theConst.tGridSize, theConst.mGridSize);
    metalArr_h = theConst.metalArr;
	tempArr_h = theConst.tempArr;
	emmArr_h = theConst.emmArr;
	energyArr_h = energyArrInit(theConst.binCenterSize);
	integral_h = new double[theConst.nPixX*theConst.nPixY*theConst.binCenterSize];
	
	////////////create temperature and metallicity grids for interpolation within the device
	tempGrid_h = (double*) calloc(theConst.tGridSize, sizeof(double));
    tempGrid_h = tempGridInit(theConst.tGridSize, theConst.tempAxis);
	metalGrid_h = (double*) calloc(theConst.mGridSize, sizeof(double));
    metalGrid_h = metalGridInit(theConst.mGridSize, theConst.metalAxis);
    //    printf("tempGrid[%d] = %f\n", theConst.tGridSize-1, tempGrid_h[theConst.tGridSize-1]);
    //    printf("metalGrid[%d] = %f\n", theConst.mGridSize-1, metalGrid_h[theConst.mGridSize -1]);
    rotMat_h = (double*) calloc(9,sizeof(double));
    rotMat_h = rotMatInit(theConst.theta, theConst.phi, theConst.epsilon);
    size_t sizeRebin = theConst.tGridSize*theConst.mGridSize*theConst.binCenterSize*sizeof(double);
    eLoadTime = clock();
    
    //create space for device tensors and transfer data from host to device tensors
	printf("Allocating memory to GPU tensors and loading into GPU\n");
    sGPUTime = clock();
    cudaMalloc((void**)&metalArr_d, sizeT);
	cudaMemcpy(metalArr_d, metalArr_h, sizeT, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&tempArr_d, sizeT);
	cudaMemcpy(tempArr_d, tempArr_h, sizeT, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&emmArr_d, sizeT);
	cudaMemcpy(emmArr_d, emmArr_h, sizeT, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&integral_d, sizeInt);
	cudaMemcpy(integral_d, integral_h, sizeInt, cudaMemcpyHostToDevice);
    cudaMalloc((void**)&rebinArr_d, sizeRebin);
	cudaMemcpy(rebinArr_d, rebinArr_h, sizeRebin, cudaMemcpyHostToDevice);
    cudaMalloc((void**)&tempGrid_d, sizeTGrid);
    cudaMemcpy(tempGrid_d, tempGrid_h, sizeTGrid, cudaMemcpyHostToDevice);
    cudaMalloc((void**)&metalGrid_d, sizeMGrid);
	cudaMemcpy(metalGrid_d, metalGrid_h, sizeMGrid, cudaMemcpyHostToDevice);
    cudaMalloc((void**)&rotMat_d, 9*sizeof(double));
	cudaMemcpy(rotMat_d, rotMat_h,  9*sizeof(double), cudaMemcpyHostToDevice);
    eGPUTime = clock();
    
    ///////////////test to see if cuda node is functioning
    double *testArr_h = new double[10];
    testArr_h[0] = 0.0;
    double *testArr_d;
    bool cudaWorking = false;
    size_t testSize = 10*sizeof(double);
    cudaMalloc((void**)&testArr_d, testSize);
	cudaMemcpy(testArr_d, testArr_h, testSize, cudaMemcpyHostToDevice);
    cudFunctionTest<<<1,10>>>(testArr_d);
    cudaMemcpy(testArr_h, testArr_d, testSize, cudaMemcpyDeviceToHost);
    if(testArr_h[0]>0.0){
        cudaWorking = true;
        printf("Cuda node is functional\n");
    }else{
        printf("Cuda node is broken. Fix before proceeding.\n");
    }
    
    bool debugging = false;
    //////////////call the integration
    printf("Running integration\n");
    sIntTime = clock();
    double * combinedRegion;
    if (cudaWorking) {
        integrate2<<<theConst.nPixX*theConst.nPixY, theConst.binCenterSize>>>(rebinArr_d, tempArr_d, metalArr_d, emmArr_d, integral_d, tempGrid_d, metalGrid_d, rotMat_d,theConst, debugging);
        printf("Transferring integration results back to CPU\n \n");
        cudaMemcpy(integral_h, integral_d, sizeInt, cudaMemcpyDeviceToHost);
        eIntTime = clock();
        
        //convert integral_h back to a tensor
        double *** integralMatrix = makeIntegralMatrix(integral_h, theConst.nPixX, theConst.nPixY, theConst.binCenterSize); 
        bool displayRegion = true;
        ///////////////////display the spectra at different combined locations
        if(displayRegion){
            int gridX1 = x1;
            int gridX2 = x2;
            int gridY1 = y1;
            int gridY2 = y2;
            combinedRegion = sumArea(gridX1,gridX2,gridY1,gridY2,integralMatrix,theConst);
//            for (int i=0; i<theConst.binCenterSize; i++) {
//                printf("bin[%d] = %f  %f  ",i, combinedRegion[i], log10(powf(10,integralMatrix[1][2][i]) + powf(10,integralMatrix[2][2][i])));
//            }
        }
        
        //////////////print the spectra
        bool printSpect = false;
        if(printSpect){
            printSpectra(integralMatrix, theConst.nPixX, theConst.nPixY);
        }
        eTime = clock();
        bool printClocking = true;
        if (printClocking) {
            printf("\nNumber of integrations: %d\n", theConst.nPixX*theConst.nPixY*theConst.binCenterSize);
            printf("Time to load input tensors: %f seconds\n", (double(eLoadTime-sLoadTime))/CLOCKS_PER_SEC);
            printf("Time to process cooling function and rebin: %f seconds\n", (double(eCoolTime-sCoolTime))/CLOCKS_PER_SEC);
            printf("Time to interate: %f seconds\n", (double(eIntTime-sIntTime))/CLOCKS_PER_SEC);
            printf("Time to allocate space and transfer tensors to GPU: %f seconds\n", (double(eGPUTime-sGPUTime))/CLOCKS_PER_SEC);
            printf("Total run time: %f seconds\n", (double(eTime-sTime))/CLOCKS_PER_SEC);
            
        }
        //////////////debugging integral function cuda call
        if(debugging){
//            integrate<<<theConst.nPixX*theConst.nPixY, theConst.binCenterSize>>>(rebinArr_d, tempArr_d, metalArr_d, emmArr_d, integral_d, tempGrid_d, metalGrid_d, theConst, debugging);
//            
//            cudaMemcpy(integral_h, integral_d, sizeInt, cudaMemcpyDeviceToHost);
            for (int i=0; i<11; i++) {
                printf("integral[%d] = %f\n", i, integral_h[i]);
            }
        }
    }
    
    ////////////testing new functions
    bool testing = false;
    if(testing){
//        for (int i=0; i<5; i++) {
//            printf("TempArr[%d][%d][%d] = %f\n", i,i,i,tenRetrieveH(tempArr_h, theConst.nx, theConst.ny, theConst.nz, i, i, i));
//            tempArr_h[i*theConst.ny*theConst.nz + i*theConst.nz + i] = 0.0;
//            printf("0000--TempArr[%d][%d][%d] = %f\n", i,i,i,tenRetrieveH(tempArr_h, theConst.nx, theConst.ny, theConst.nz, i, i, i));
//            
//        }
//        tempInit<<<theConst.nx*theConst.ny, theConst.nz>>>(function2, tempArr_d, theConst);
//        cudaMemcpy(tempArr_h, tempArr_d, sizeT, cudaMemcpyDeviceToHost);
//        for (int i=0; i<5; i++) {
//            printf("TempArr[%d][%d][%d] = %f\n", i,i,i,tenRetrieveH(tempArr_h, theConst.nx, theConst.ny, theConst.nz, i, i, i));
//        }
        bool blocking = false;
        
        //test for blocking
        if (blocking) {
            for (int i=1; i<2; i++) {
                sleep(2);
                printf("integration #%d\n", i);
                integrate<<<theConst.nPixX*theConst.nPixY, theConst.binCenterSize>>>(rebinArr_d, tempArr_d, metalArr_d, emmArr_d, integral_d, tempGrid_d, metalGrid_d, theConst, debugging);
            }
            cudaMemcpy(integral_h, integral_d, sizeInt, cudaMemcpyDeviceToHost);
        }
    }    
    return combinedRegion;
	
}

int main(){
    jaco_state ja;
    ja.nlastbin = 256;
    ja.rshock = 1.0;
    ja.angdist = 10;
    ja.pixscale = 10;
    ja.nx = 128;
    ja.ny = 128;
    ja.theta = 0;
    ja.phi = 0;
    ja.epsilon = 0;
    ja.a_ell = 1.0;
    ja.b_ell = 1.0;
    ja.n_ell = 256*256*256;
    
	file_reader coolingFile;
	double centBin[ja.nlastbin];
	strcpy(coolingFile.spectralCode, "/home/tchap/mekal.bin");
	coolingFile.debug = 2;
	read_cooling_function(coolingFile);
    makeReBin(centBin, ja.nlastbin);
	rebincoolingfunction(centBin, ja.nlastbin, coolingFile);  
    //double float issue
    ja.rebinnedcooling = coolingFile.rebinnedCooling; 
    ja.tgridsize = coolingFile.tGridSize;
    ja.mgridsize = coolingFile.mGridSize;
    ja.egridsize = coolingFile.eGridSize;
    ja.tempaxis = coolingFile.tempAxis;
    ja.metalaxis = coolingFile.metalAxis;
    ja.Tprof = tempArrInit(ja.n_ell, ja.a_ell, ja.b_ell, ja.rshock);
    ja.Zprof = metalArrInit(ja.n_ell, ja.a_ell, ja.b_ell, ja.rshock);
    ja.nenhprof = emmArrInit(ja.n_ell, ja.a_ell, ja.b_ell, ja.rshock);
    double *results = runSimulation(ja, 2,4,3,5);
    for (int i=0; i<ja.nlastbin; i++) {
        printf("bin[%d] = %f  ",i, results[i]);
    }
    Cleanup();
	return 0;
}




