#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "config.h"
#include "debug.h"

void cpyBack(){
	HANDLE_ERROR(cudaMemcpy(hVel,d_hVel,sizeof(vector3) * NUMENTITIES, cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(hPos,d_hPos,sizeof(vector3) * NUMENTITIES, cudaMemcpyDeviceToHost));
}

__global__ void computeAccel(vector3 ** d_accels, vector3 * d_hPos, double * d_mass){
	int j,k;
	int i = blockIdx.x*blockDim.x+threadIdx.x;
	if(i>=NUMENTITIES)return;
	for (j=0;j<NUMENTITIES;j++){
		if (i==j) {
			FILL_VECTOR(d_accels[i][j],0,0,0);
		}
		else{
			vector3 distance;
			for (k=0;k<3;k++) distance[k]=d_hPos[i][k]-d_hPos[j][k];
			double magnitude_sq=distance[0]*distance[0]+distance[1]*distance[1]+distance[2]*distance[2];
			double magnitude=sqrt(magnitude_sq);
			double accelmag=-1*GRAV_CONSTANT*d_mass[j]/magnitude_sq;
			FILL_VECTOR(d_accels[i][j],accelmag*distance[0]/magnitude,accelmag*distance[1]/magnitude,accelmag*distance[2]/magnitude);
		}
	}
}

__global__ void sumAccel(vector3 * d_hVel, vector3 * d_hPos, vector3 ** d_accels){
	int j,k;
	int i = blockIdx.x*blockDim.x+threadIdx.x;
	if(i>=NUMENTITIES)return;

	//copied from second for loop in compute
	vector3 accel_sum={0,0,0};

	for (j=0;j<NUMENTITIES;j++){
		for (k=0;k<3;k++)
			accel_sum[k]+=d_accels[i][j][k];
	}
	//compute the new velocity based on the acceleration and time interval
	//compute the new position based on the velocity and time interval
	for (k=0;k<3;k++){
		d_hVel[i][k]+=accel_sum[k]*INTERVAL;
		d_hPos[i][k]+=d_hVel[i][k]*INTERVAL;
	}
}

//compute: Updates the positions and locations of the objects in the system based on gravity.
//Parameters: None
//Returns: None
//Side Effect: Modifies the hPos and hVel arrays with the new positions and accelerations after 1 INTERVAL
void compute(){
	//make an acceleration matrix which is NUMENTITIES squared in size;

	//vector3* values=(vector3*)malloc(sizeof(vector3)*NUMENTITIES*NUMENTITIES);
	//values is a 1d array, accels is a way to access values with 2d syntax
	//vector3** accels=(vector3**)malloc(sizeof(vector3*)*NUMENTITIES);
	// for (i=0;i<NUMENTITIES;i++)
	// 	accels[i]=&values[i*NUMENTITIES];

	//first compute the pairwise accelerations.  Effect is on the first argument.
	//for loop for computing accelerations, which are then stored in the values array
	//IN THEORY: start a kernel with i threads and compute and sum up the elements of the row there
	computeAccel<<<NUMENTITIES/BLOCKSIZE+1,BLOCKSIZE>>>(d_accels, d_hPos, d_mass);

	/*for (i=0;i<NUMENTITIES;i++){
		//kernel needed on this outer for loop

		for (j=0;j<NUMENTITIES;j++){
			if (i==j) {
				FILL_VECTOR(accels[i][j],0,0,0);
			}
			else{
				vector3 distance;
				for (k=0;k<3;k++) distance[k]=hPos[i][k]-hPos[j][k];
				double magnitude_sq=distance[0]*distance[0]+distance[1]*distance[1]+distance[2]*distance[2];
				double magnitude=sqrt(magnitude_sq);
				double accelmag=-1*GRAV_CONSTANT*mass[j]/magnitude_sq;
				FILL_VECTOR(accels[i][j],accelmag*distance[0]/magnitude,accelmag*distance[1]/magnitude,accelmag*distance[2]/magnitude);
			}
		}
	}*/
	//sum up the rows of our matrix to get effect on each entity, then update velocity and position.
	sumAccel<<<NUMENTITIES/BLOCKSIZE+1,BLOCKSIZE>>>(d_hVel, d_hPos, d_accels);
	//for (i=0;i<NUMENTITIES;i++){
		//kernel needed on this outer for loop
		// vector3 accel_sum={0,0,0};

		// for (j=0;j<NUMENTITIES;j++){
		// 	for (k=0;k<3;k++)
		// 		accel_sum[k]+=accels[i][j][k];
		// }
		// //compute the new velocity based on the acceleration and time interval
		// //compute the new position based on the velocity and time interval
		// for (k=0;k<3;k++){
		// 	hVel[i][k]+=accel_sum[k]*INTERVAL;
		// 	hPos[i][k]+=hVel[i][k]*INTERVAL;
		// }
	//}
	// free(accels);
	// free(values);
	cpyBack();
}
