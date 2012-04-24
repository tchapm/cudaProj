//1d Integration for an N number of barriers

#include <iostream>
#include <math.h>
#include <time.h>
#include <stdio.h>
using namespace std;
int functCall=0;

struct pixel{
	float temp;
	float metal;
	float parts;
	//float lamda;
	float integral;
	int E;
};

struct grid{
	int x;
	int y;
};

float f(int x, float k){
	float value;
	//value = 2*x;
	value= float(powf(powf(1/(k+powf(x,2)), 3.0)+powf(1/(k+powf(x,2)), 5.0/2.0), 0.5));
	functCall++;
	return value;
}


void varCalc(pixel* list, grid* xyCoords, int listSize, int gridSize, float err){
	float k;
	float h, tFunct;
	float nextVal, actErr;
	int a, b, step=1;
	a=0;
	b=10;
	float last;
	tFunct=0.5*(f(a,k)+f(b,k));
	double prevStep[200];
	prevStep[0]=tFunct;
	int n=1;
	int jump= listSize/gridSize;
	
	for (int j=0; j<gridSize; j++) {
		for (int i=jump*j; i<(jump*(j+1)); i++)  {  //section for f(x)
			k=powf(xyCoords[j].x,2)+powf(xyCoords[j].y,2);
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
}

void init(grid inputData[], pixel eVal[], int nGrid, int nPixels)
{
	srand(1203);
	int i, j, k;
	k=0;
	
	while(k<nGrid){
		for(i = 10; i >= -10; i--) 
		{
			for(j=10; j>=-10; j--){
				inputData[k].y = j;
				inputData[k].x = i;
				k++;
			}
		}
	}
	j=1;
	for (i=1; i < nPixels; i++) {
		if (j>256) {
			j=1;
		}
		eVal[i].E = j;
		j++;
	}
	
}



int main (int argc, char* argv[]) {	
	float accuracy=0.0001;
	int numPixels=512*512;
	int numGrid=21*21;
	int i,j;
	clock_t sTime, eTime;
	pixel pixelArr[numPixels];
	grid gridArr[numGrid];
	
	
	
	init(gridArr, pixelArr, numGrid, numPixels);
	sTime = clock();
	varCalc(pixelArr, gridArr, numPixels, numGrid, accuracy);
	eTime = clock();
	int step= numPixels/numGrid;
	
	/*for (j=0; j<numGrid; j++) {
	 for (i=(step*j); i<(step*(j+1)); i++) {
	 printf("For pixel: %d ", i);
	 printf("X: %d Y: %d E: %d ", gridArr[j].x, gridArr[j].y, pixelArr[i].E);
	 printf("Integral: %f \n\n", pixelArr[i].integral);
	 
	 }
	 }*/
	printf("Total run time: %f seconds\n", (double(eTime-sTime))/CLOCKS_PER_SEC);
	printf("Number of integrations: %d\n", numPixels);
	
	return 0;
}



