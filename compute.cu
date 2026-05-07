/*
Parallelised by Samhain Ackerman and Henry Leap
*/

#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "config.h"

//compute: Updates the positions and locations of the objects in the system based on gravity.
//Parameters: None
//Returns: None
//Side Effect: Modifies the hPos and hVel arrays with the new positions and accelerations after 1 INTERVAL
__global__ void compute(
        vector3 * d_hVel,
        vector3 * d_hPos,
        double * d_mass,
        vector3 ** d_accels
){
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
