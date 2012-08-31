
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

#define PI 3.14159265

using namespace std;

double *metalArr_h;
double *metalArr_d;
double *emmArr_h;
double *emmArr_d;
double *tempArr_h;
double *tempArr_d;
double *energyArr_h;
double *ellArr_h;
double *ellArr_d;
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
constants theConst;


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
    theConst.psi = ja.psi;
    theConst.a_ell = ja.a_ell;
    theConst.b_ell = ja.b_ell;
    theConst.n_ell = ja.n_ell;
//    theConst.rMax = ja.rshock;
    theConst.rebinnedCooling = ja.rebinnedcooling;
    theConst.metalArr=ja.Zprof;
    theConst.tempArr=ja.Tprof;
    theConst.emmArr=ja.nenhprof;
	theConst.accuracy=0.0001;
	theConst.nz = (int)ja.rshock; //this should be changed!!! to theConst.rMax?
    theConst.ellMax = ja.ell_max;//powf(theConst.nPixX*theConst.nPixX + theConst.nPixY*theConst.nPixY + theConst.rMax*theConst.rMax,0.5); //????
    
}
//jaco_state ja;
//ja.nlastbin = 256;
//ja.rshock = 2.0;
//ja.angdist = 10;
//ja.pixscale = 10;
//ja.nx = 128;
//ja.ny = 128;
//ja.theta = 0;
//ja.phi = 0;
//ja.epsilon = 0;
//ja.a_ell = 1.0;
//ja.b_ell = 1.0;
//ja.n_ell = ja.nx*ja.ny*ja.rshock;

double*** runSimulation(jaco_state ja){
    ////////////set constants and sizes of data structures
    clock_t sTime, eTime, sLoadTime, eLoadTime, sGPUTime, eGPUTime, sCoolTime, eCoolTime, sIntTime, eIntTime;
    /////////////start clock
    sTime = clock();
	
	constInit(theConst, ja);
    
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
//    ellArr_h = ellArrInit(theConst.n_ell,2*ja.nx*ja.ny);
    
	integral_h = new double[theConst.nPixX*theConst.nPixY*theConst.binCenterSize];
	
	////////////create temperature and metallicity grids for interpolation within the device
	tempGrid_h = (double*) calloc(theConst.tGridSize, sizeof(double));
    tempGrid_h = tempGridInit(theConst.tGridSize, theConst.tempAxis);
	metalGrid_h = (double*) calloc(theConst.mGridSize, sizeof(double));
    metalGrid_h = metalGridInit(theConst.mGridSize, theConst.metalAxis);
    //    printf("tempGrid[%d] = %f\n", theConst.tGridSize-1, tempGrid_h[theConst.tGridSize-1]);
    //    printf("metalGrid[%d] = %f\n", theConst.mGridSize-1, metalGrid_h[theConst.mGridSize -1]);
    rotMat_h = (double*) calloc(9,sizeof(double));
    rotMat_h = rotMatInit(theConst.theta, theConst.phi, theConst.psi);
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
	cudaMemcpy(rotMat_d, rotMat_h, 9*sizeof(double), cudaMemcpyHostToDevice);
//    cudaMalloc((void**)&ellArr_d, theConst.n_ell*sizeof(double));
//	cudaMemcpy(ellArr_d, ellArr_h, theConst.n_ell*sizeof(double), cudaMemcpyHostToDevice);
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
    double*** integralMatrix;
    if (cudaWorking) {
        integrate2<<<theConst.nPixX*theConst.nPixY, theConst.binCenterSize>>>(rebinArr_d, tempArr_d, metalArr_d, emmArr_d, integral_d, tempGrid_d, metalGrid_d, rotMat_d, theConst, debugging);
        printf("Transferring integration results back to CPU\n \n");
        cudaMemcpy(integral_h, integral_d, sizeInt, cudaMemcpyDeviceToHost);
        eIntTime = clock();
        
        //convert integral_h back to a tensor
        integralMatrix = makeIntegralMatrix(integral_h, theConst.nPixX, theConst.nPixY, theConst.binCenterSize); 
        
        //////////////print the spectra
        bool printSpect = false;
        if(printSpect){
            printSpectra(integralMatrix, energyArr_h, theConst);
        }
        bool sumSpect = false;
        if(sumSpect){
//            plotImage(integralMatrix, energyArr_h, theConst);
//            sumSpectra(integralMatrix, energyArr_h, theConst);
        }
        eTime = clock();
        bool printClocking = false;
        if (printClocking) {
            printf("\nNumber of integrations: %d\n", theConst.nPixX*theConst.nPixY*theConst.binCenterSize);
            printf("Time to load input tensors: %f seconds\n", (double(eLoadTime-sLoadTime))/CLOCKS_PER_SEC);
            printf("Time to process cooling function and rebin: %f seconds\n", (double(eCoolTime-sCoolTime))/CLOCKS_PER_SEC);
            printf("Time to integrate: %f seconds\n", (double(eIntTime-sIntTime))/CLOCKS_PER_SEC);
            printf("Time to allocate space and transfer tensors to GPU: %f seconds\n", (double(eGPUTime-sGPUTime))/CLOCKS_PER_SEC);
            printf("Total run time: %f seconds\n", (double(eTime-sTime))/CLOCKS_PER_SEC);
            
        }
        //////////////debugging integral function cuda call
        if(debugging){
//            integrate<<<theConst.nPixX*theConst.nPixY, theConst.binCenterSize>>>(rebinArr_d, tempArr_d, metalArr_d, emmArr_d, integral_d, tempGrid_d, metalGrid_d, theConst, debugging);
//            
//            cudaMemcpy(integral_h, integral_d, sizeInt, cudaMemcpyDeviceToHost);
            for (int i=0; i<14; i++) {
                printf("integral[%d] = %f\n", i, integral_h[i]);
                printf("integral[%d] = %f\n", i+2080512, integral_h[i+2080512]);
            }
//            for (int i=2080512; i<2080526; i++) {
//                printf("integral[%d] = %f\n", i, integral_h[i]);
//            }
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
    return integralMatrix;
	
}

double *sumArea(int x1, int x2, int y1, int y2, double*** completeArr, int nBins){
    double *combinedRegion;
    combinedRegion = (double*) calloc(nBins, sizeof(double));
    for (int k=0; k<nBins; k++) {
        for (int i=x1; i<=x2; i++) {
            for (int j=y1; j<=y2; j++) {
                combinedRegion[k]+= powf(10,completeArr[i][j][k]);
            }
        }
        combinedRegion[k]=log10(combinedRegion[k]);
    }
    return combinedRegion;
}

double **multiSpec(double ***totalSpecta,int inputRegions[][4], int numSpectra, int numBins){
    double **tempTotal = sci_fmatrixD(numSpectra, numBins);
    double *oneRunResults; 
    for (int i=0; i<numSpectra; i++) {
        oneRunResults = sumArea(inputRegions[i][0], inputRegions[i][1], inputRegions[i][2], inputRegions[i][3], totalSpecta, numBins);
        for(int j=0; j<numBins; j++){
            tempTotal[i][j] = oneRunResults[j]; 
        }
    }
    return tempTotal;
}

int main(){
    double degVal=90;
    jaco_state ja;
    ja.nlastbin = 256;
    ja.rshock = 2.0;
    ja.angdist = 10;
    ja.pixscale = 10;
    ja.nx = 128;
    ja.ny = 128;
    ja.phi = 0;//degVal*PI/180;
    ja.theta = degVal*PI/180;
    ja.psi = degVal*PI/180;
    ja.a_ell = 10.0;
    ja.b_ell = 1.0;
    ja.n_ell = ja.nx; //or 256*256*256?
    ja.ell_max = ja.rshock + 0.1;
    
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
    ja.Tprof = tempArrInit(ja.n_ell, ja.a_ell, ja.b_ell, ja.ell_max); //is nPixX correct?
    ja.Zprof = metalArrInit(ja.n_ell, ja.a_ell, ja.b_ell, ja.ell_max);
    ja.nenhprof = emmArrInit(ja.n_ell, ja.a_ell, ja.b_ell, ja.ell_max);
    
    //input region values as a matrix with the desired num of spectra as first component
    //the second component is the four corners of the rectangular region
    int numSpectra;
    numSpectra = 3;
    //example of input with three regions
    int inputRegions[3][4] = 
    {
        {3, 3, 3, 3, }, // row 0
        {20, 21, 20, 21, }, // row 1
        {51, 51, 51, 51} // row 2
    };
//    int anArray[3][5] =
//    {
//        { 1, 2, 3, 4, 5, }, // row 0
//        { 6, 7, 8, 9, 10, }, // row 1
//        { 11, 12, 13, 14, 15 } // row 2
//    };
    //example: to get spectra for x=10, y=15 loop through i=0->i=nlastbin for totalSpectra[10][15][i]
//    double ***totalSpectra = runSimulation(ja); //spectra for all the regions
//    char* filePrefix = "/home/tchap/NVIDIA_GPU_COMPUTING_SDK/C/src/cud3Dsim/fluxValues/";
    for(int i=10;i<=180;i+=10){
        ja.psi = i*PI/180;
        double ***totalSpectra = runSimulation(ja);
        char filePrefix[265] = "/home/tchap/NVIDIA_GPU_COMPUTING_SDK/C/src/cud3Dsim/fluxValues/";
        ostringstream convert; 
        convert << i;       
        char* inputName = new char [convert.str().size()+1];
        strcpy(inputName,convert.str().c_str());
        strcat(filePrefix,inputName);
        strcat(filePrefix,"epsilon.txt");
       
        printf("printing to %s.\n",filePrefix);
        plotImage(totalSpectra, energyArr_h, theConst, filePrefix);
    }
    //the expected output of my program as a double** where the first component is the index of the desired region and the second is the array of the spectra
    /*
    double ***totalSpectra = runSimulation(ja);
    double **collectedResults = multiSpec(totalSpectra, inputRegions, numSpectra, ja.nlastbin);
    
    double *results = sumArea(3, 3, 3, 3, totalSpectra, ja.nlastbin);
    //
    for (int i=0; i<ja.nlastbin; i++) {
//        printf("area1[%d] = %f  ",i, results[i]);
//        printf("area2[%d] = %f  ",i, collectedResults[1][i]);
//        printf("area3[%d] = %f  ",i, collectedResults[2][i]);
    }
    */
    Cleanup();
	return 0;
}




