/*
Parallelised by Samhain Ackerman and Henry Leap
*/

#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "config.h"

__device__ void reduce(vector3 accel_sum) {
	int inBlockI = threadIdx.x;
	int iSection = threadIdx.y;
	__shared__ vector3 theseRows[IS_PER_BLOCK][THREADS_PER_I];
        COPY_VECTOR(theseRows[inBlockI][iSection], accel_sum);
        FILL_VECTOR(accel_sum, 0,0,0);
        __syncthreads();
        // ignore that we only need to do this for iSection == 0
        if (iSection) return; // or don't
        for (int j = 0; j < THREADS_PER_I; j++)
                ADD_VECTORS(accel_sum, theseRows[inBlockI][j]);
}

//compute: Updates the positions and locations of the objects in the system based on gravity.
//Parameters: None
//Returns: None
//Side Effect: Modifies the hPos and hVel arrays with the new positions and accelerations after 1 INTERVAL
__global__ void computeVel(
        vector3 * d_hVel,
        vector3 * d_hPos,
        double * d_mass
){
	int j,k;
	int i = blockIdx.x*blockDim.x+threadIdx.x;
	int incr = blockDim.y, offset = threadIdx.y;
	vector3 accel_sum={0,0,0};
	if(i>=NUMENTITIES) goto reduce;
	//if(i>=NUMENTITIES) goto pos;


	for (j = offset; j < NUMENTITIES; j += incr){
	        vector3 accel;
		if (i==j) {
			FILL_VECTOR(accel,0,0,0);
		}
		else{
			vector3 distance;
			for (k=0;k<3;k++) distance[k]=d_hPos[i][k]-d_hPos[j][k];
			double magnitude_sq=distance[0]*distance[0]+distance[1]*distance[1]+distance[2]*distance[2];
			double magnitude=sqrt(magnitude_sq);
			double accelmag=-1*GRAV_CONSTANT*d_mass[j]/magnitude_sq;
			FILL_VECTOR(accel,accelmag*distance[0]/magnitude,accelmag*distance[1]/magnitude,accelmag*distance[2]/magnitude);
		}
		for (k=0;k<3;k++) accel_sum[k]+=accel[k];
	}

	reduce:

	reduce(accel_sum);

	if (i >= NUMENTITIES || offset) return;
	// __shared__ vector3 theseRows; // would be nice to declare down here
	// COPY_VECTOR(theseRows[blockIdx.x][offset], accel_sum);
	// for (int l = SAME_I_THREADS >> 1; l; l >>= 1) {
	//        __syncthreads();
	//        if (offset >= l) continue;
	//        ADD_VECTORS(
	//                theseRows[blockIdx.x][offset],
	//                theseRows[blockIdx.x][offset + l]
	//        );
	//}
        // __syncthreads();

	// if(i >= NUMENTITIES || offset != 0) return;

        // for (int l = 1; l < SAME_I_THREADS; l++)
        //        ADD_VECTORS(accel_sum, theseRows[blockIdx.x][l]);

	// COPY_VECTOR(accel_sum, theseRows[blockIdx.x][0]);

	//compute the new velocity based on the acceleration and time interval
	//compute the new position based on the velocity and time interval
	//for (k=0;k<3;k++){
		// d_hVel[i][k] =
	//	atomicAdd(&d_hVel[i][k], accel_sum[k]*INTERVAL);
		// d_hPos[i][k]+=d_hVel[i][k]*INTERVAL;
	//}
	//pos:
	//__syncthreads();
	//if(i >= NUMENTITIES || offset != 0) return;
	for (k=0;k<3;k++)
		d_hVel[i][k]+=accel_sum[k]*INTERVAL;

}

__global__ void updatePos(vector3 * d_hVel, vector3 * d_hPos) {
	int k, i = blockIdx.x*blockDim.x+threadIdx.x;

	if (i >= NUMENTITIES) return;

	for (k=0;k<3;k++)
		d_hPos[i][k]+=d_hVel[i][k]*INTERVAL;
}

void compute(vector3 * d_hVel, vector3 * d_hPos, double * d_mass) {
	dim3 threadsPerBlock(IS_PER_BLOCK, THREADS_PER_I);
        computeVel<<<GRIDSIZE,threadsPerBlock>>>(d_hVel, d_hPos, d_mass);
        cudaDeviceSynchronize();
        updatePos<<<GRIDSIZE,IS_PER_BLOCK>>>(d_hVel, d_hPos);
}
