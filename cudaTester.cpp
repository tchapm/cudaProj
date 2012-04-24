#include <iostream>
#include <vector>
#include <stdlib.h>
#include <fcntl.h>
#include <errno.h>
#include <math.h>
#include <sstream>
//#include <sciutils.h>
//#include <jaco.h>
#include "data_structs.h"

using namespace std;

struct file_reader{
	int debug;
	char spectralCode[51];
	//string spectralCode;
	int eGridSize;
	int tGridSize;
	int mGridSize;
	float *** masterCooling;
	float *** rebinnedCooling;
	float ** tempMetalMat;
	
	vector<double> energyAxis;
	vector<double> tempAxis;
	vector<double> metalAxis;
	vector<double> lastbin;
	
	/*double* energyAxis;
     double* tempAxis;
     double* metalAxis;
     double* lastbin;
     */
	float* rebinArr;
	double opz;
	double eBinSize;
	long nLastBin;
	int l1, l2;
	
};
struct pixel{
	float temp;
	float metal;
	float emm;
	float ** tempMat;
	float ** metalMat;
	float ** emmMat;
	float *** tempTens;
	float *** metalTens;
	float *** emmTens;
	//float lamda;
	float integral;
	int E;
};

struct constants {
	float depth; // in Mpc
	int nPix;
	float pixToMpc;
	int binCenterSize;
	float accuracy;
	int nGrid;
	int nx;
	int ny;
	int nz;
	int eGridSize;
	int tGridSize;
	int mGridSize;
	
};

float* metalGrid;
float* tempGrid;


void read_cooling_function(struct file_reader& js)
{
	long k,j;
	int i;
	
	int coolfile;
	
	//printf("Taking H_0 = %d km/s/Mpc, Flat Universe, Omega_matter = %.2f\n",
	// HUBBLE,OMEGA_M);
	
	if (js.debug > 1)
		cout << "Reading cooling function from " << js.spectralCode << endl;
	
	// Now read in the description of the cooling function.
	if ((coolfile = open(js.spectralCode,O_RDONLY)) < 0){
		cout << "Could not open cooling file " << js.spectralCode << endl;
		perror("errno");
		exit(1);
	}
	
	// The first row contains the grid size in energy, temperature, and
	// metallicity.
	
	
	read(coolfile,&js.eGridSize,sizeof(js.eGridSize));
	read(coolfile,&js.tGridSize,sizeof(js.tGridSize));
	read(coolfile,&js.mGridSize,sizeof(js.mGridSize));
	//printf("%d %d %d\n",js.eGridSize,js.tGridSize,js.mGridSize);
	
	// Allocate memory to contain the whole cooling function
	//js->mastercooling = sci_ftensor(js->tGridSize,js->mGridSize,js->eGridSize);
	/*js.masterCooling = new double *[js.tGridSize];
	 for (i=0; i<js.tGridSize; i++) {
	 js.masterCooling[i] = new double[js.mGridSize];
	 }
	 */
	js.masterCooling = sci_ftensor(js.tGridSize,js.mGridSize,js.eGridSize);
	
	//js.energyAxis = (double*) malloc(js.eGridSize);
	//js.tempAxis = (double*) malloc(js.tGridSize);
	//js.metalAxis = (double*) malloc(js.mGridSize);
	js.energyAxis.reserve(js.eGridSize);
	js.tempAxis.reserve(js.tGridSize);
	js.metalAxis.reserve(js.mGridSize);
	
	// Next we read in the centers of the bins over which the cooling
	// function was integrated. This appears only once.
	
	
	for(i=0; i<js.eGridSize; i++)
		read(coolfile, &js.energyAxis[i],sizeof(double)); 
	for(i=0; i<js.tGridSize; i++)
		read(coolfile, &js.tempAxis[i],sizeof(double));
	for(i=0; i<js.mGridSize; i++)
		read(coolfile, &js.metalAxis[i],sizeof(double));
	
	js.eBinSize = js.energyAxis[2]-js.energyAxis[1]; 	
	// Now we read the temperature and metallicity abscissae, followed
	// by the cooling function itself. Note that the former repeat
	// themselves redundantly throughout the coolfile; preceding each
	// array (of length egridsize) containing the cooling function,
	// there are two numbers specifying the temperature and metallicity
	// abscissae assigned to that particular function, and these must
	// self-consistently yield a uniform grid in (T,Z) space.
	
	int l;
	for (k = 0; k < js.tGridSize; k++){
		for (j = 0; j < js.mGridSize; j++){
			for (l=0; l< js.eGridSize; l++){
				read(coolfile,&js.masterCooling[k][j][l],sizeof(float));
			}
		}
	}			   
	
	
	close(coolfile);
	
	if (js.debug > 0) {
		printf("Done. Cooling function has %d x %d x %d elements\n",
			   js.tGridSize,js.mGridSize,js.eGridSize);
		
		printf("Energy range: %E to %E keV\n", 
			   js.energyAxis[0], js.energyAxis[js.eGridSize-1]);  
		
		printf("Temperature range: %E to %E keV\n",
			   pow(10, js.tempAxis[0]), pow(10, js.tempAxis[js.tGridSize-1]));
		
		printf("Metallicity range: %E to %E keV\n",
			   pow(10, js.metalAxis[0]), pow(10, js.metalAxis[js.mGridSize-1]));
		
		printf("Cooling function range: %E to %E\n", 
			   js.masterCooling[0][0][0], js.masterCooling[js.tGridSize-1][js.mGridSize-1][js.eGridSize-1]);
		
		
	}
	
	//js->coolaccel = gsl_interp_accel_alloc(); ???
}

void printCooling(struct file_reader& js){
	/*for (k=(js.tGridSize-1)/4; k<js.tGridSize; k+=(js.tGridSize-1)/4) {
	 printf("tempVal: %E \t metalVal: %E \n", js.tempAxis[k], js.metalAxis[k]);
	 }*/
	
	printf("ListPlot[{{");
	int k;
	for (k=0; k<js.eGridSize-1; k++) {
		printf("{%f, %f},", js.energyAxis[k], log10(js.masterCooling[(js.tGridSize-1)/4][(js.mGridSize-1)/4][k]));
	}
	printf("{%f, %f}},{", js.energyAxis[k], log10(js.masterCooling[(js.tGridSize-1)/4][(js.mGridSize-1)/4][k]));
	
	for (k=0; k<js.eGridSize-1; k++) {
		printf("{%f, %f},", js.energyAxis[k], log10(js.masterCooling[(js.tGridSize-1)*2/4][(js.mGridSize-1)*2/4][k]));
	}
	printf("{%f, %f}},{", js.energyAxis[k], log10(js.masterCooling[(js.tGridSize-1)*2/4][(js.mGridSize-1)*2/4][k]));
	
	for (k=0; k<js.eGridSize-1; k++) {
		printf("{%f, %f},", js.energyAxis[k], log10(js.masterCooling[(js.tGridSize-1)*3/4][(js.mGridSize-1)*3/4][k]));
	}
	printf("{%f, %f}}}]", js.energyAxis[k], log10(js.masterCooling[(js.tGridSize-1)*3/4][(js.mGridSize-1)*3/4][k]));
	
	/*for (k=0; k<js.eGridSize; k++) {
	 printf("%.3E\n", js.energyAxis[k]);
	 }
	 
	 for (k = 3; k < 4; k++){
	 for (j = 0; j < 50; j++){
	 printf("%.3E ",js.masterCooling[k][j][js.eGridSize/4]);
	 printf("\n");
	 }
	 
	 }/*
	 printf("\n\n");
	 for (k = js.tGridSize*4/5; k < js.tGridSize; k++){
	 for (j = js.mGridSize*4/5; j < js.mGridSize; j++){
	 printf("%.3E ",js.masterCooling[k][j][js.eGridSize-1]);
	 }
	 printf("\n");
	 }*/
	
	/*for (i=0; i<js.eGridSize/2; i++) {
	 printf("%.2E\t %.2E\n", js.masterCooling[1][1][i], js.masterCooling[2][2][i]);
	 }*/
	
	/*double ave=0;
	 for (l=0; l<js.eGridSize; l++) {
	 for(k=0; k<js.tGridSize; k++){
	 for (j=0; j<js.mGridSize; j++) {
	 ave+=js.masterCooling[k][j][l];
	 }
	 }
	 //printf("%.0E %.0E   ",js.masterCooling[24][24][l],js.energyAxis[l]);
	 printf("%.6E", ave/(js.mGridSize*js.tGridSize));
	 ave=0;
	 printf("\n");
	 }*/
}


void printRebin(struct file_reader& js, double *binCenter){
	int k;
	/*printf("ListPlot[{{");
     for (k=0; k<js.eGridSize-1; k++) {
     printf("{%f, %f},", log10(js.energyAxis[k]), log10(js.masterCooling[(js.tGridSize-1)*3/4][(js.mGridSize-1)*3/4][k]));
     }
     printf("{%f, %f}},{", log10(js.energyAxis[k]), log10(js.masterCooling[(js.tGridSize-1)*3/4][(js.mGridSize-1)*3/4][k]));
     for (k=0; k<255; k++) {
     printf("{%f, %f},", log10(binCenter[k]), (js.rebinnedCooling[k][(js.tGridSize-1)*3/4][(js.mGridSize-1)*3/4]));
     }
     printf("{%f, %f}}}]", log10(binCenter[k]), (js.rebinnedCooling[k][(js.tGridSize-1)*3/4][(js.mGridSize-1)*3/4]));
     */
	printf("ListPlot[{{");
	for (k=0; k<255; k++) {
		printf("{%f, %f},", log10(binCenter[k]), js.rebinnedCooling[k][32][40]);
	}
	printf("{%f, %f}},{", log10(binCenter[k]), js.rebinnedCooling[k][32][40]);
	for (k=0; k<255; k++) {
		printf("{%f, %f},", log10(binCenter[k]), tenRetrieveH(js.rebinArr, 256, js.tGridSize, js.mGridSize, k, 32, 40));
	}
	printf("{%f, %f}}}]", log10(binCenter[k]), tenRetrieveH(js.rebinArr, 256, js.tGridSize, js.mGridSize, k, 32, 40));
    
	/*printf("\n\ncooling function = ");
	 for (k=0; k<10; k++) {
	 printf("{%f, %f}\n", (js.energyAxis[k]), (js.masterCooling[36][36][k]));
	 }
	 printf("rebinned = ");
	 for (k=0; k<5; k++) {
	 printf("{%f, %f}\n", (binCenter[k]), (js.rebinnedCooling[k][36][36]));
	 }*/
	
	/*printf("ListPlot[{");
     for (k=0; k<255; k++) {
     printf("{%f, %f},", log10(binCenter[k]), js.rebinnedCooling[k][32][40]);
     }
     printf("{%f, %f}}]", log10(binCenter[k]), js.rebinnedCooling[k][32][40]);
     */
	/*printf("ListPlot[{");
	 for (k=0; k<255; k++) {
	 if (js.rebinnedCooling[k][36][36]< -1) {
	 cout << endl << js.rebinnedCooling[k][36][36] << "\t k = " << k;
	 //printf("{%f, %f},", log10(binCenter[k]), (js.rebinnedCooling[k][(js.tGridSize-1)*3/4][(js.mGridSize-1)*3/4]));
	 }
	 }*/
	//printf("{%f, %f}}]", log10(binCenter[k]), (js.rebinnedCooling[k][(js.tGridSize-1)*3/4][(js.mGridSize-1)*4/4]));
}


void makeReBin(double *binCenter, int numElem){
	//is binCenter the center of each rebinned bin???
	double binSize = 12.0/numElem;
	binCenter[0] = 0.1;
	for (int k = 1; k<256; k++) {
		binCenter[k] = binSize*k + binCenter[0];
	}
}

int binSearch(struct file_reader& js, double e){
	int mid, min = 0;
	int max = js.eGridSize-1; 
	double mval, mplus;
	int pos;
	while(min < max){
		mid = (min+max)/ 2;
		mval = js.energyAxis[mid];
		if (e >= mval){
			mplus = js.energyAxis[mid+1];
			if (e < mplus){
				pos = mid;
				break;
			}
			min = mid + 1;
		}
		else
			max = mid;
	}
	//if (min < max)
	//	pos = js.energyAxis[js.eGridSize-1];
	return pos;
}

// Take the finely-grained cooling function in mastercooling
// and rebin it to match the detector binning given in binlo and binhi
int rebincoolingfunction(double *binCenter, int nBin, struct file_reader& js)
{
	long      i,j,k;
	unsigned long e1pos,e2pos;
	double    mincool = 1.E30,e1fac,e2fac,e1,e2,logmincool=-1E30,countBinSize;
	int      haschanged=1,l;
	
	js.eBinSize = js.energyAxis[1]-js.energyAxis[0]; //I added and am unsure???
	
	//js->l1 = js->l2 = -1;
	if (haschanged) {
		// Generate a new coolingfunction matrix
		js.l1 = 0; js.l2 = nBin-1;
		js.rebinnedCooling = sci_ftensor(nBin,js.tGridSize,js.mGridSize);
		js.lastbin.reserve(nBin*sizeof(double));
		//js.lastbin = (double*) malloc(nBin*sizeof(double));
		printf("Rebinning cooling function: %d bins requested.\n",nBin);
		
		for  (l = 0; l < nBin; ++l)
			js.lastbin[l] = binCenter[l];
		js.nLastBin = nBin;
	} else
		// We have checked that rebinnedcooling[1..nbin] must be a subset
		// of the previously calculated rebinnedcooling[1..nlastbin];
		// specifically, it equals rebinnedcooling[l1..l2]. Time to exit.
		return 0;
	
	// We assume that all the bin sizes are the SAME, and that
	// the grid is REGULAR.
	countBinSize = (binCenter[1]-binCenter[0]);
	if (js.debug > 0)
		printf("spectral bin size is %f keV\n",countBinSize);
	if (js.debug > 1)
		printf("Bins start at %f and end at %f\n", binCenter[0],binCenter[nBin-1]);
	
	//cout << "\ntest cooling j = 49 k=48 e2pos= 146" << js.masterCooling[49][48][146];
	for (l = 0; l < nBin; ++l) {
		// Calculate the rest-frame bin corresponding to the given
		// energy range
		js.opz = 1.0; //correction for z
		e1 = js.opz*(binCenter[l]);
		e2 = js.opz*(binCenter[l]+countBinSize);
		
		
		if (e1 < js.energyAxis[0] || e2 < js.energyAxis[0] ||
			e1 > js.energyAxis[js.eGridSize-1] ||
			e2 > js.energyAxis[js.eGridSize-1]) {
			printf("For bounds: %f--%f\n",e1,e2);
			printf("Energy bounds of cooling matrix were exceeded.\n");
		}
		
		// Find the position of the bins corresponding to this range
		/*e1pos = gsl_interp_accel_find(js.coolaccel,
		 js.energyaxis,js.egridsize,e1);
		 e2pos = gsl_interp_accel_find(js.coolaccel,
		 js.energyaxis,js.egridsize,e2);
		 */
		//bisection or faster search for later	
		e1pos = binSearch(js, e1);
		e2pos = binSearch(js, e2);
		/*for (k=0; k<(js.eGridSize-1); k++) {
		 //if (e1>=js.energyAxis[k] && e1<js.energyAxis[k+1]) {
		 //	e1pos=k; //or maybe js.energyAxis[k]
		 //}
		 if (e2>=js.energyAxis[k] && e2<js.energyAxis[k+1]) {
		 e2pos=k; //or maybe js.energyAxis[k]
		 }
		 }*/
		//cout << "\ne1 = "<< e1 << "\te2 = " << e2 << endl;
		if (e1pos > e2pos){
			printf("Lower energy bound %ld > upper energy bound %ld.", e1pos , e2pos);
			exit(-1);
		}
		
		if (e1pos != e2pos) {
			// Upper and lower energy limits occur in different bins
			e1fac = (js.energyAxis[e1pos+1]-e1);
			e2fac = (e2-js.energyAxis[e2pos]);
		} else {
			// Upper and lower energy limits occur in the same bin
			// (not desireable, since this means the cooling function is not
			// binned with fine enough resolution)
			e1fac = (e2-e1);
			e2fac = 0;
		}
		
		// Add together the contributions from the starting and ending bins
		for (j = 0; j < js.tGridSize; ++j){
			for (k = 0; k < js.mGridSize; ++k) {
				js.rebinnedCooling[l][j][k] =
				js.masterCooling[j][k][e1pos]+ //took out *e1fac
				js.masterCooling[j][k][e2pos]*e2fac;
				/*if (l==14 && j == 35 && k == 35){
				 cout << "\nthe value of rebincooling = " << log10(js.rebinnedCooling[l][j][k]) << endl;
				 cout << "e1pos = " << e1pos << "\t e2pos = " << e2pos << "\t e1fac = "
				 << e1fac << "\te2fac = " << e2fac << endl;
				 }
				 if (l==1 && j == 36 && k == 36){
				 cout << "\n1. the value of rebincooling = " << (js.rebinnedCooling[l][j][k]) << endl;
				 cout << "e1pos = " << e1pos << "\t e2pos = " << e2pos << "\t e1fac = "
				 << e1fac << "\te2fac = " << e2fac << "   e1 = " << e1 << "   e2  = " << e2 << endl;
				 cout << "Master cooling " << js.masterCooling[j][k][e1pos] << "\t" << js.masterCooling[j][k][e2pos] << endl;
				 }*/
				if (js.rebinnedCooling[l][j][k] > 0. &&
					js.rebinnedCooling[l][j][k] < mincool)
					mincool = js.rebinnedCooling[l][j][k];
			}
		}
		
		// Add together the contributions from the intermediate bins
		for (i = (long)e1pos+1; i < (long)e2pos; ++i){
			for (j = 0; j < js.tGridSize; ++j)
				for (k = 0; k < js.mGridSize; ++k) {
					js.rebinnedCooling[l][j][k] +=
					js.masterCooling[j][k][i]*js.eBinSize;
					/*if (l==1 && j == 36 && k == 36){
					 cout << "\n2. the value of rebincooling = " << (js.rebinnedCooling[l][j][k]) << endl;
					 cout << "e1pos = " << e1pos << "\t e2pos = " << e2pos << "\t e1fac = "
					 << e1fac << "\te2fac = " << e2fac << endl;
					 cout << "master cooling " << log10(js.masterCooling[j][k][i]) << "\t BinSize " << js.eBinSize << endl;
					 }*/
					if (js.rebinnedCooling[l][j][k] > 0. &&
						js.rebinnedCooling[l][j][k] < mincool)
						mincool = js.rebinnedCooling[l][j][k];
				}
		}
		
		logmincool = log10(mincool);
		
		// Now take the log of the rebinned function, for interpolation in
		// log space. Truly null cooling rates are simply set to the
		// minimum rate seen altogether.
		for (j = 0; j < js.tGridSize; ++j){
			for (k = 0; k < js.mGridSize; ++k) {
				js.rebinnedCooling[l][j][k] =
				(js.rebinnedCooling[l][j][k] < mincool ? logmincool :
				 log10(js.rebinnedCooling[l][j][k]));
			}
		}
		
	}
	return 1;
}
void inputSet(float **inputData, int eBin, int tGridSize, int mGridSize, float ***masterCooling){
	int i, j;
	for(i = 0; i < tGridSize; i++) 
	{
		for(j=0; j< mGridSize; j++){
			inputData[i][j] = log10(masterCooling[i][j][eBin]);
		}
	}
	
}
/*void tempMatSet(struct file_reader& js, pixel eVal[]){
 int i, j, eBin=5;
 //eVal[eBin].tempMetalMat=sci_fmatrix(js.tGridSize, js.mGridSize);
 for(eBin = 0; eBin<js.eGridSize; eBin++){ 
 
 eVal[eBin].tempMetalMat=sci_fmatrix(js.tGridSize, js.mGridSize);
 for(i = 0; i < js.tGridSize; i++) 
 {
 for(j=0; j< js.mGridSize; j++){
 
 eVal[eBin].tempMetalMat[i][j] = log10(js.masterCooling[i][j][eBin]);
 }
 }
 }
 }*/

void sphereInit(pixel& init){
	int x, y, z;
	init.tempTens=sci_ftensor(256, 256, 256);
	init.metalTens=sci_ftensor(256, 256, 256);
	init.emmTens=sci_ftensor(256, 256, 256);
	float r;
    
	for (x=0; x<256; x++) {
		for(y=0; y<256; y++){
			for (z=0; z<256; z++) {
				r=powf((powf(x-128,2)+powf(y-128,2)+powf(z-128, 2)), -0.5);
				
				init.metalTens[x][y][z] = powf((r+1), -0.5);
				init.tempTens[x][y][z] = 6*init.metalTens[x][y][z];
				init.emmTens[x][y][z] = powf((pow(r,2)+1), -3);//emmission measure
			}
		}
	}
	//cout << init.tempMat[0][0][0];
}

/*float bilinInterp(float* interpMat, float x, float y, constants theConst, int blockId){
 int x1, x2, y1, y2;
 double temp1, temp2, R1, R2, interpValue;
 double ten11, ten12, ten21, ten22;
 x1 = (int)x;
 x2 = x1+1;
 y1 = (int)y;
 y2 = y1+1;
 int nx = theConst.nx;
 int ny = theConst.ny;
 int nz = theConst.binCenterSize;
 
 temp1 = (x2-x)/(x2-x1);
 temp2 = (x-x1)/(x2-x1);
 ten11 = get3DValFrom1DArray(interpMat, nx, ny, nz, x1, y1, blockId);
 ten21 = get3DValFrom1DArray(interpMat, nx, ny, nz, x2, y1, blockId);
 ten12 = get3DValFrom1DArray(interpMat, nx, ny, nz, x1, y2, blockId);
 ten22 = get3DValFrom1DArray(interpMat, nx, ny, nz, x2, y2, blockId);
 
 ten11 = tenRetrieveH(interpMat, nx, ny, nz, x1, y1, blockId);
 ten21 = tenRetrieveH(interpMat, nx, ny, nz, x2, y1, blockId);
 ten12 = tenRetrieveH(interpMat, nx, ny, nz, x1, y2, blockId);
 ten22 = tenRetrieveH(interpMat, nx, ny, nz, x2, y2, blockId);
 //printf("temp1: %f temp2: %f", temp1, temp2);
 //printf("coord(1,1), (1,2), (2,1), (2,2): %f %f %f %f\n", interpMat[x1*ny + y1], interpMat[x1*ny + y2], interpMat[x2*ny + y1], interpMat[x2*ny + y2]);
 //R1 = temp1*interpMat[x1*ny + y1] + temp2*interpMat[x2*ny + y1];
 //R2 = temp1*interpMat[x1*ny + y2] + temp2*interpMat[x2*ny + y2];
 /*R1 = temp1*ten11 + temp2*ten21;
 R2 = temp1*ten12 + temp2*ten22;
 interpValue = ((y2-y)/(y2-y1))*R1 + ((y-y1)/(y2-y1))*R2;
 
 return interpValue;
 }*/

int getLowerIndex(float* axisArr, float value, int arrLength){
	for(int i=0; i<arrLength; i++){
		if (value<axisArr[i]) {
			return i-1;
		}
	}
	return arrLength-2;
}


/*
 float f(int x, float k){
 float value;
 //value = 2*x;
 value= float(powf(powf(1/(k+powf(x,2)), 3.0)+powf(1/(k+powf(x,2)), 5.0/2.0), 0.5));
 functCall++;
 return value;
 }
 */

float tenRetrieveD(float* oneDArray, int nx, int ny, int nz, int x, int y, int z){
	int index = x*ny*nz + y*nz + z;
	return oneDArray[index];
}
float bilinInterpVal(float* interpMat, float y, float z, int energyBin, float* tempGrid, float* metalGrid, constants theConst){
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
	ten11 = powf(10,tenRetrieveD(interpMat, nx, ny, nz, energyBin, y1, z1));
	ten21 = powf(10,tenRetrieveD(interpMat, nx, ny, nz, energyBin, y2, z1));
	ten12 = powf(10,tenRetrieveD(interpMat, nx, ny, nz, energyBin, y1, z2));
	ten22 = powf(10,tenRetrieveD(interpMat, nx, ny, nz, energyBin, y2, z2));
	R1 = temp1*ten11 + temp2*ten21;
	R2 = temp1*ten12 + temp2*ten22;
	interpValue = ((z2val-z)/(z2val-z1val))*R1 + ((z-z1val)/(z2val-z1val))*R2;
	
	return interpValue;
}

float linearInterpZ(float* interpMat, int x, int y, float z, constants theConst){
	int z1, z2;
	int nx = theConst.nx;
	int ny = theConst.ny;
	int nz = theConst.nz;
	double ten1, ten2, interpValue;
	z1 = (int)z;
	z2 = z1+1;
	ten1 = tenRetrieveD(interpMat, nx, ny, nz, x, y, z1);
	ten2 = tenRetrieveD(interpMat, nx, ny, nz, x, y, z2);
	interpValue = ten1 + (ten2 -ten1)*(z - z1);
	return interpValue;
	
}

float integrate(float* rebinCool, float* tempArr, float* metalArr, float* emmArr, float* tempGrid, float* metalGrid, constants theConst, int threadIndex, int blockIndex){
	//integrate from depth/2 to -depth/2
    //	int i=blockDim.x * blockIdx.x + threadIdx.x; //threadIdx.x is the channel/energy-bin	
	int j= blockIndex;   
	int x, y;
	x = j/theConst.nx;
	y = j%theConst.nx;
	int energyBin = threadIndex;
	int a, b, n, step=1;
	float tFunct, h, actErr, last, nextVal;
	float T, Z, rebinA, rebinB;
	double prevStep[200];
	
	//b = theConst.depth/2 + theConst.nz/2;
	//a = theConst.nz/2 - theConst.depth/2;
	b = theConst.depth + theConst.nz/2;
	a = theConst.nz/2;
	h = b-a;
	actErr = 1.0;
	last = 1.E30;
	nextVal = 0.0;
	n = 1;
	T = tenRetrieveD(tempArr, theConst.nx, theConst.ny, theConst.nz, x, y, b);
	Z = tenRetrieveD(metalArr, theConst.nx, theConst.ny, theConst.nz, x, y, b);
    
	rebinB = bilinInterpVal(rebinCool, T, Z, energyBin, tempGrid, metalGrid, theConst);
	T = tenRetrieveD(tempArr, theConst.nx, theConst.ny, theConst.nz, x, y, a);
	Z = tenRetrieveD(metalArr, theConst.nx, theConst.ny, theConst.nz, x, y, a);
    
	rebinA = bilinInterpVal(rebinCool, T, Z, energyBin, tempGrid, metalGrid, theConst);
	tFunct = 0.5*powf(tenRetrieveD(emmArr, theConst.nx, theConst.ny, theConst.nz, x, y, a), 2.0)*rebinA;
	tFunct += 0.5*powf(tenRetrieveD(emmArr, theConst.nx, theConst.ny, theConst.nz, x, y, b), 2.0)*rebinB;
	//tFunct += 0.5*rebinB;
    
    
	while (actErr>=0.1) {
		step=step*2;
		h= float(b-a)/step;		
		for (int l=1; l<step; l=l+2) {
			T = linearInterpZ(tempArr, x, y, l*h+a, theConst);
			Z = linearInterpZ(metalArr, x, y, l*h+a, theConst);
			nextVal+=bilinInterpVal(rebinCool, T, Z, energyBin, tempGrid, metalGrid, theConst);
		}
		nextVal+=prevStep[n-1];
		prevStep[n]=nextVal;
		nextVal=h*(nextVal);
		actErr=fabs(last-nextVal);
		last=nextVal;
        //        printf("LastVal = %f\n",last);
		n++;
	}
	return last;
	//TODO change to actual indices to make more readable for Andi in future
	//integral[i] = last;	 
}
/*void integrate(float* rebinCool, float* tempArr, float* metalArr, float* emmArr, float* energyArr, float* integral, constants theConst){
 //integrate from depth/2 to -depth/2
 
 //when do i need to call the bilinear interpolation? 
 int x, y, z;
 
 int a, b, n, step=1;
 float tFunct, h, actErr, last, nextVal;
 float T, Z, rebinA, rebinB;
 double prevStep[200];
 
 b = theConst.depth/2;
 a = -b;
 h = b-a;
 actErr = 1.0;
 last = 1.E30;
 nextVal = 0.0;
 n = 1;
 
 for (int j=0; j<gridSize; j++) {
 x = j/theConst.nx;
 y = j%theConst.nx;
 for (int i=jump*j; i<(jump*(j+1)); i++){
 T = get3DValFrom1DArray(tempArr, theConst.nx, theConst.ny, theConst.nz, x, y, a+128);
 Z = get3DValFrom1DArray(metalArr, theConst.nx, theConst.ny, theConst.nz, x, y, a+128);
 rebinA = bilinInterp(rebinCool, T, Z, theConst, j);
 T = get3DValFrom1DArray(tempArr, theConst.nx, theConst.ny, theConst.nz, x, y, b+128);
 Z = get3DValFrom1DArray(metalArr, theConst.nx, theConst.ny, theConst.nz, x, y, b+128);
 rebinB = bilinInterp(rebinCool, T, Z, theConst, j);
 
 tFunct = 0.5*powf(get3DValFrom1DArray(emmArr, theConst.nx, theConst.ny, theConst.nz, x, y, a+128), 2.0)*rebinA;
 tFunct += 0.5*powf(get3DValFrom1DArray(emmArr, theConst.nx, theConst.ny, theConst.nz, x, y, b+128), 2.0)*rebinB;
 prevStep[0]=tFunct;
 
 while (actErr>theConst.accuracy) {
 step=step*2;
 h= float(b-a)/step;
 
 for (int l=1; l<step; l=l+2) {
 T = get3DValFrom1DArray(tempArr, theConst.nx, theConst.ny, theConst.nz, x, y, l*h+a+128);
 Z = get3DValFrom1DArray(metalArr, theConst.nx, theConst.ny, theConst.nz, x, y, l*h+a+128);
 nextVal+=bilinInterp(rebinCool, T, Z, theConst, j);
 }
 nextVal+=prevStep[n-1];
 prevStep[n]=nextVal;
 nextVal=h*(nextVal);
 actErr=fabs(last-nextVal);
 last=nextVal;
 n++;
 }	
 
 /*float *rebinMatrix = new float[theConst.nx*theConst.ny];
 for (int m=0; m<nx; m++) {
 for (int n=0; n<ny; n++) {
 rebinMatrix[m*ny + n] = rebinCool[m*ny*nz + n*nz + i];
 }
 } 
 //x*ny*nz + y*nz + i for index of rebinCool
 tFunct = 0.5;*//*
                 //coolingFunction = bilinInterp(rebinCool, T, Z, theConst.ny);
                 //TODO change to actual indices to make more readable for Andi in future
                 integral[i] = tFunct;
                 }
                 }
                 
                 }*/
/*void varCalc(pixel* list, struct file_reader& js, int listSize, float err){
 float k;
 float h, tFunct;
 float nextVal, actErr;
 
 int a, b, step=1;
 a=0;
 b=10;
 float last;
 //tFunct=0.5*(f(a,k)+f(b,k));
 double prevStep[200];
 prevStep[0]=tFunct;
 int n=1;
 int gridSize = js.eGridSize; //??? 
 int jump= listSize/gridSize;
 
 for (int j=0; j<gridSize; j++) {
 for (int i=jump*j; i<(jump*(j+1)); i++)  {  //section for f(x)
 //use energy levels of [j] and [j+1]
 
 //k=powf(xyCoords[j].x,2)+powf(xyCoords[j].y,2);
 h=(b-a);
 actErr=1.0;
 last=1.E30;
 nextVal=0.0;
 n=1;
 tFunct=0.5*f(a, k);
 tFunct+=0.5*f(b, k); 
 prevStep[0]=tFunct;
 while (actErr>=err) { //will break when difference of
 //successive calculations less than err
 step=step*2;
 h= float(b-a)/step;
 nextVal=0.0;
 for (int l=1; l<step; l=l+2) {
 nextVal+=f(l*h+a, k);
 }
 nextVal+=prevStep[n-1];
 prevStep[n]=nextVal;
 nextVal=h*(nextVal);
 actErr=fabs(last-nextVal);
 //cout << "last " << last << endl;
 last=nextVal;
 n++;
 }		
 //cout << "n: " << n << endl;
 //cout << "last: " << last << endl;
 list[i].integral= float(last*powf(list[i].E,-0.5));
 //cout <<  i << endl << list[i].E << endl;
 //cout << list[i].integral << endl;
 prevStep[0]=tFunct;
 n=1;
 step=1;
 }
 }	
 }*/

void constInit(constants &theConst){
	theConst.depth = 2.0; //2.0  in Mpc
	theConst.nPix = 128*128;
	theConst.pixToMpc = 0.01;
	theConst.binCenterSize = 256;
	theConst.accuracy=0.0001;
	theConst.nGrid=1;
	theConst.nx=256;
	theConst.ny=256;
	theConst.nz=256;
}
int main (int argc, char * const argv[]) {
    file_reader coolingFile;
	constants theConst;
	constInit(theConst);
	int rebinSize = theConst.binCenterSize;
	int nx, ny, nz;
	nx=10; ny=20; nz=30;
	float ***flatTen = sci_ftensor(nx, ny, nz);
	float *tenArr = new float[nx*ny*nz];
	double centBin[rebinSize]; 
	flatTen[2][10][13] = 5.4;
	tenArr = tensorTo1DArray(flatTen, nx, ny, nz);
	tenArr[621] = 2.2;
	//float *rebinArr;
	float *metalArr = metalArrInit(nx, ny, nz);
	float *emmArr = emmArrInit(nx, ny, nz);
	float *tempArr = tempArrInit(nx, ny, nz);
	float *energyArr = energyArrInit(theConst.binCenterSize);
	float *integralArr = new float[theConst.nPix*theConst.binCenterSize];
    //	
    //	printf("metalArr = %f\n", get3DValFrom1DArray(metalArr, 1, 2, 3, nx, ny, nz));
    //	printf("emmArr = %f\n", get3DValFrom1DArray(emmArr, 1, 2, 3, nx, ny, nz));
    //	printf("tempArr = %f\n", get3DValFrom1DArray(tempArr, 1, 2, 3, nx, ny, nz));
    
	//printf("flatTen[1][2][3] = %f\t tenArr = %f\n", flatTen[2][10][13], get3DValFrom1DArray(tenArr, 2, 10, 13, nx, ny, nz));
	//int k;
	makeReBin(centBin, theConst.binCenterSize);
	strcpy(coolingFile.spectralCode, "/Users/tchap/Desktop/cudaResearch/mekal.bin");
	coolingFile.debug=2;
	read_cooling_function(coolingFile);
	theConst.eGridSize = coolingFile.eGridSize;
	theConst.tGridSize = coolingFile.tGridSize;
	theConst.mGridSize = coolingFile.mGridSize;
    //	printf("const:   binCenterSize = %d \ttGridSize = %d \t mGridSize = %d", theConst.binCenterSize, theConst.tGridSize, theConst.mGridSize);
    //	printf("\ncoolFile:   binCenterSize = %d \ttGridSize = %d \t mGridSize = %d", theConst.binCenterSize, coolingFile.tGridSize, coolingFile.mGridSize);
	
	rebincoolingfunction(centBin, theConst.binCenterSize, coolingFile);
	coolingFile.rebinArr = tensorTo1DArray(coolingFile.rebinnedCooling, theConst.binCenterSize, coolingFile.tGridSize, coolingFile.mGridSize);
	tempGrid = (float*) calloc(theConst.tGridSize, sizeof(float));
	for(int i =0; i<theConst.tGridSize; i++){
		tempGrid[i] = pow(10,coolingFile.tempAxis[i]);
		//printf("\ntempGrid[%d] = %f", i, tempGrid[i]);
	}
	metalGrid = (float*) calloc(theConst.mGridSize, sizeof(float));
	for (int i = 0; i<theConst.mGridSize; i++) {
		metalGrid[i] = pow(10,coolingFile.metalAxis[i]);
		//printf("\nmetalGrid[%d] = %f", i, metalGrid[i]);
	}
	
	//bilinInterp(rebinArr, 2.5, 0.3, theConst, 1);
    //	printf("\nrebinArr[150][2][3] = %f", log10(bilinInterpCool(coolingFile.rebinArr, 3.0, 1.0, theConst, 150)));
    //	printf("\nrebinArr[150][2][3] = %f", coolingFile.rebinnedCooling[150][31][39]);
    //	printf("\nrebin1D[150][2][3] = %f\n", tenRetrieveH(coolingFile.rebinArr, rebinSize, coolingFile.tGridSize, coolingFile.mGridSize, 150, 32, 40));
	
	/*for (int i = 0; i<256; i++){
     coolingFile.rebinnedCooling[i][][];
     }*/
	
	//bilinInterpCool(rebinArr, 2.5, 0.3, theConst, 1);
	//printf("rebinArr[16382][2.5][0.3] = %f", bilinInterpCool(rebinArr, 2.5, 0.3, theConst, 163));
	//printf("interpolation[0][0][0]= %f\tactual[0][0][0]= %f", bilinInterp(rebinArr, 0.0, 0.0, 0), h_pixelArr[0].tempMat[0][0][0]);
    
	//float ** inputData = sci_fmatrix(newFile.tGridSize, newFile.mGridSize);
	//inputSet(inputData, k, newFile.tGridSize, newFile.mGridSize, newFile.masterCooling);
	//cout << "test: " << inputData[0][0];
	//cout << "\nmasterCool: " << log10(newFile.masterCooling[0][0][1]);
	//printCooling(newFile);
	/*for(int i =0; i<50; i++){
     printf("\ntemp[%d] = %f\n", i, pow(10,coolingFile.tempAxis[i]));
     printf("metal[%d] = %f\n", i, pow(10,coolingFile.metalAxis[i]));
     }*/
	//printRebin(coolingFile, centBin);
	//int numPixels=newFile.eGridSize;
	//cout << newFile.masterCooling[0][0][0];
	//pixel h_pixelArr[numPixels];//=(pixel*)calloc(numPixels, sizeof(pixel));
	//tempMatSet(newFile, h_pixelArr);
	//sphereInit(h_pixelArr[0]);
	//cout << h_pixelArr[0].tempMat[0][0][0];
	//printf("interpolation[1.2][1.9][0]= %f\n", bilinInterp(tempArr, 1.2, 1.9, ny));
	//printf("interpolation[0][0][0]= %f\tactual[0][0][0]= %f", bilinTempInterpolation(h_pixelArr[0], 0.0, 0.0, 0), h_pixelArr[0].tempMat[0][0][0]);
    //theConst.nPix, theConst.binCenterSize
    int threadIndex = 1;
    int blockIndex = 1;
    float integrateVal = integrate(coolingFile.rebinArr, tempArr, metalArr, emmArr, tempGrid, metalGrid, theConst, threadIndex, blockIndex);
    printf("IntegrationVal = %f\n", integrateVal);
    return 0;
}
