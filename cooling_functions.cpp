#include "data_structs.cpp"

using namespace std;

struct file_reader{
	int debug;
	char spectralCode[51];
	//string spectralCode;
	int eGridSize;
	int tGridSize;
	int mGridSize;
	float *** masterCooling;
	float *** rebinnedCooling; //possibly should be double
	float * rebinArr;
//	vector<double> energyAxis;
//	vector<double> tempAxis;
//	vector<double> metalAxis;
//	vector<double> lastbin;
    double *energyAxis;
	double *tempAxis;
	double *metalAxis;
	double *lastbin;
	double opz;
	double eBinSize;
	long nLastBin;
	int l1, l2;
	
};

void read_cooling_function(struct file_reader& js) {
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
	
	
//	js.energyAxis.reserve(js.eGridSize);
//	js.tempAxis.reserve(js.tGridSize);
//	js.metalAxis.reserve(js.mGridSize);
    js.energyAxis = new double[js.eGridSize];
    js.tempAxis = new double[js.tGridSize];
    js.metalAxis = new double[js.mGridSize];
	
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
	
//	if (js.debug > 0) {
		printf("Done. Cooling function has %d x %d x %d elements\n",
			   js.tGridSize,js.mGridSize,js.eGridSize);
		
		printf("Energy range: %E to %E keV\n", 
			   js.energyAxis[0], js.energyAxis[js.eGridSize-1]);  
		
		printf("Temperature range: %E to %E keV\n",
			   js.tempAxis[0],js.tempAxis[js.tGridSize-1]);
		
		printf("Metallicity range: %E to %E keV\n",
			   js.metalAxis[0],js.metalAxis[js.mGridSize-1]);
//		
//		printf("Cooling function range: %E to %E\n", 
//			   js.masterCooling[0][0][0], js.masterCooling[js.tGridSize-1][js.mGridSize-1][js.eGridSize-1]);
//		
//		
//	}
	
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
	 
	 }
	 printf("\n\n");
	 for (k = js.tGridSize*4/5; k < js.tGridSize; k++){
	 for (j = js.mGridSize*4/5; j < js.mGridSize; j++){
	 printf("%.3E ",js.masterCooling[k][j][js.eGridSize-1]);
	 }
	 printf("\n");
	 }
	
     for (i=0; i<js.eGridSize/2; i++) {
	 printf("%.2E\t %.2E\n", js.masterCooling[1][1][i], js.masterCooling[2][2][i]);
	 }
	
	double ave=0;
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
	 for (k=0; k<10; k++) {e
	 printf("{%f, %f}\n", (js.energyAxis[k]), (js.masterCooling[36][36][k]));
	 }
	 printf("rebinned = ");
	 for (k=0; k<5; k++) {
	 printf("{%f, %f}\n", (binCenter[k]), (js.rebinnedCooling[k][36][36]));
	 }*/
	
	/*printf("ListPlot[{");
	 for (k=0; k<255; k++) {
	 printf("{%f, %f},", log10(binCenter[k]), (js.rebinnedCooling[k][(js.tGridSize-1)*3/4][(js.mGridSize-1)*3/4]));
	 }
	 printf("{%f, %f}}]", log10(binCenter[k]), (js.rebinnedCooling[k][(js.tGridSize-1)*3/4][(js.mGridSize-1)*4/4]));
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
	float binSize = 12.0/numElem;
	binCenter[0] = 0.1;
    int k;
	for (k = 1; k<255; k++) {
		binCenter[k] = binSize + binCenter[k-1];
//		printf("%f, ", binCenter[k]);
	}
    binCenter[k] = binSize + binCenter[k-1];

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
        js.lastbin = new double[nBin];
//		js.lastbin.reserve(nBin*sizeof(double));
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
//	if (js.debug > 0)
//		printf("spectral bin size is %f keV\n",countBinSize);
//	if (js.debug > 1)
//		printf("Bins start at %f and end at %f\n", binCenter[0],binCenter[nBin-1]);
	
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
				js.masterCooling[j][k][e1pos]*e1fac+ //took out *e1fac but then replaced it
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



