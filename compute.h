/*
Parallelised by Samhain Ackerman and Henry Leap
*/

__global__ void compute(
        vector3 * d_hVel,
        vector3 * d_hPos,
        double * d_mass,
        vector3 ** d_accels
);
