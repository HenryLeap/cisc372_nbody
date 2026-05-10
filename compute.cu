/*
Parallelised by Samhain Ackerman and Henry Leap
*/

#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "config.h"

__device__ void reduce(vector3 accel_sum) {
        // NOTE: This `offset` is *very* different from the one in `computeVel`
        for (int offset = WARPSIZE/2; offset; offset >>= 1) {
                accel_sum[0] += __shfl_down_sync(0xFFFFFFFF, accel_sum[0], offset);
                accel_sum[1] += __shfl_down_sync(0xFFFFFFFF, accel_sum[1], offset);
                accel_sum[2] += __shfl_down_sync(0xFFFFFFFF, accel_sum[2], offset);
        }
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
	int i = blockIdx.x;
	int offset = threadIdx.x;
	vector3 accel_sum={0,0,0};

	for (j = offset; j < NUMENTITIES; j += BLOCKSIZE){
	        vector3 accel;
		if (i==j) continue;
                vector3 displacement;
                for (k = 0; k < 3; k++)
                        displacement[k] = d_hPos[i][k] - d_hPos[j][k];
                double distance_sq = (
                        displacement[0] * displacement[0] +
                        displacement[1] * displacement[1] +
                        displacement[2] * displacement[2]);
                double distance = sqrt(distance_sq);
                double accelmag = -1 * GRAV_CONSTANT * d_mass[j] / distance_sq;
                FILL_VECTOR(accel,
                        accelmag * displacement[0] / distance,
                        accelmag * displacement[1] / distance,
                        accelmag * displacement[2] / distance);
		ADD_VECTORS(accel_sum, accel);
	}


	reduce(accel_sum);

	if (offset) return;

	ADD_VECTORS(d_hVel[i], INTERVAL * accel_sum);
}

__global__ void updatePos(vector3 * d_hVel, vector3 * d_hPos) {
	int i = blockIdx.x*blockDim.x+threadIdx.x;

	if (i >= NUMENTITIES) return;

	ADD_VECTORS(d_hPos[i], INTERVAL * d_hVel[i]);
}

void compute(vector3 * d_hVel, vector3 * d_hPos, double * d_mass) {
        computeVel<<<GRIDSIZE,BLOCKSIZE>>>(d_hVel, d_hPos, d_mass);
        cudaDeviceSynchronize();
        updatePos<<<GRIDSIZE,BLOCKSIZE>>>(d_hVel, d_hPos);
}
