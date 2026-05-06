#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "vector.h"
#include "config.h"
#include "planets.h"
#include "compute.h"

// represents the objects in the system.  Global variables
vector3 *hVel, *d_hVel;
vector3 *hPos, *d_hPos;
vector3 **d_accels;
vector3 *d_values;
double *mass, *d_mass;

//initHostMemory: Create storage for numObjects entities in our system
//Parameters: numObjects: number of objects to allocate
//Returns: None
//Side Effects: Allocates memory in the hVel, hPos, and mass global variables
void initHostMemory(int numObjects)
{
	hVel = (vector3 *)malloc(sizeof(vector3) * numObjects);
	hPos = (vector3 *)malloc(sizeof(vector3) * numObjects);
	mass = (double *)malloc(sizeof(double) * numObjects);
}

__global__ void initAccels(vector3 ** d_accels, vector3 * d_values){
	int i = blockIdx.x*blockDim.x+threadIdx.x;
	if (i>=NUMENTITIES) return;
	d_accels[i]=&d_values[i*NUMENTITIES];
}
void initDevMemory(int numObjects)
{
	cudaMalloc(&d_hVel, sizeof(vector3) * numObjects);
	cudaMalloc(&d_hPos, sizeof(vector3) * numObjects);
	cudaMalloc(&d_mass, sizeof(double) * numObjects);
	cudaMalloc(&d_accels, sizeof(vector3*) * numObjects);
	cudaMalloc(&d_values, sizeof(vector3) * numObjects * numObjects);
	initAccels<<<numObjects/BLOCKSIZE+1,BLOCKSIZE>>>(d_accels, d_values);
}


//freeHostMemory: Free storage allocated by a previous call to initHostMemory
//Parameters: None
//Returns: None
//Side Effects: Frees the memory allocated to global variables hVel, hPos, and mass.
void freeHostMemory()
{
	free(hVel);
	free(hPos);
	free(mass);
}

void freeDevMemory()
{
	cudaFree(d_hVel);
	cudaFree(d_hPos);
	cudaFree(d_mass);
	cudaFree(d_accels);
	cudaFree(d_values);
}

//planetFill: Fill the first NUMPLANETS+1 entries of the entity arrays with an estimation
//				of our solar system (Sun+NUMPLANETS)
//Parameters: None
//Returns: None
//Fills the first 8 entries of our system with an estimation of the sun plus our 8 planets.
void planetFill(){
	int i,j;
	double data[][7]={SUN,MERCURY,VENUS,EARTH,MARS,JUPITER,SATURN,URANUS,NEPTUNE};
	for (i=0;i<=NUMPLANETS;i++){
		for (j=0;j<3;j++){
			hPos[i][j]=data[i][j];
			hVel[i][j]=data[i][j+3];
		}
		mass[i]=data[i][6];
	}
}

//randomFill: FIll the rest of the objects in the system randomly starting at some entry in the list
//Parameters: 	start: The index of the first open entry in our system (after planetFill).
//				count: The number of random objects to put into our system
//Returns: None
//Side Effects: Fills count entries in our system starting at index start (0 based)
void randomFill(int start, int count)
{
	int i, j, c = start;
	for (i = start; i < start + count; i++)
	{
		for (j = 0; j < 3; j++)
		{
			hVel[i][j] = (double)rand() / RAND_MAX * MAX_DISTANCE * 2 - MAX_DISTANCE;
			hPos[i][j] = (double)rand() / RAND_MAX * MAX_VELOCITY * 2 - MAX_VELOCITY;
			mass[i] = (double)rand() / RAND_MAX * MAX_MASS;
		}
	}
}

//printSystem: Prints out the entire system to the supplied file
//Parameters: 	handle: A handle to an open file with write access to prnt the data to
//Returns: 		none
//Side Effects: Modifies the file handle by writing to it.
void printSystem(FILE* handle){
	int i,j;
	for (i=0;i<NUMENTITIES;i++){
		fprintf(handle,"pos=(");
		for (j=0;j<3;j++){
			fprintf(handle,"%lf,",hPos[i][j]);
		}
		printf("),v=(");
		for (j=0;j<3;j++){
			fprintf(handle,"%lf,",hVel[i][j]);
		}
		fprintf(handle,"),m=%lf\n",mass[i]);
	}
}

void cpyOver(int numObjects){
	cudaMemcpy(d_hVel,hVel,sizeof(vector3) * numObjects, cudaMemcpyHostToDevice);
	cudaMemcpy(d_hPos,hPos,sizeof(vector3) * numObjects, cudaMemcpyHostToDevice);
	cudaMemcpy(d_mass,mass,sizeof(double) * numObjects, cudaMemcpyHostToDevice);
}

int main(int argc, char **argv)
{
	clock_t t0=clock();
	int t_now;
	//srand(time(NULL));
	srand(1234);
	initHostMemory(NUMENTITIES);
	planetFill();
	initDevMemory(NUMENTITIES);
	randomFill(NUMPLANETS + 1, NUMASTEROIDS);
	cpyOver(NUMENTITIES);

	//now we have a system.
	#ifdef DEBUG
	printSystem(stdout);
	#endif
	for (t_now=0;t_now<DURATION;t_now+=INTERVAL){
		compute();
	}
	clock_t t1=clock()-t0;
#ifdef DEBUG
	printSystem(stdout);
#endif
	printf("This took a total time of %f seconds\n",(double)t1/CLOCKS_PER_SEC);

	freeHostMemory();
	freeDevMemory();
}
