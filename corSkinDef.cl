/*
	corSkinDef kernel
*/

__kernel void corSkinDef(
	__global float* finalPos,							//float3
	__global const float* initialPos,					//float3
	__global const float* weights,
	__global const float4* matrices,
	__global const float4* quaternions,
	__global const float4* CoRs,
	const uint positionCount,
	const uint numTransforms)
{
	// this is the CUDA equavalent for indexing the arrays
	unsigned int positionId = get_global_id(0);				// access finalPos and initialPos using this value
	if (positionId >= positionCount ) return;					// We create an execute unit for more indices then we have data for, just exit early if this guy if one of the extras
	
	unsigned int positionOffset = positionId * 3;				// Base positions are float3 when they come in here!
	float4 initialPosition;
	initialPosition.x = initialPos[positionOffset];
	initialPosition.y = initialPos[positionOffset+1];
	initialPosition.z = initialPos[positionOffset+2];
	initialPosition.w = 1.0f;

	float4 finalPosition;
	finalPosition.x = 0.0f;
	finalPosition.y = 0.0f;
	finalPosition.z = 0.0f;
	finalPosition.w = 1.0f;

	// for each vertex
	// work out the weighted sum of products for the quaternions components
	unsigned int weights_offset = numTransforms;
	unsigned int weights_index = positionId * weights_offset;
	float4 q1, q2;
	unsigned int quaternion_offset = 4;
	// init the first quat
	q1.x = 0.0f;
	q1.y = 0.0f;
	q1.z = 0.0f;
	q1.w = 0.0f;

	for (unsigned int j = 0; j < numTransforms; j++){
		// init the quat add operand
		q2.x = quaternions[j*quaternion_offset];
		q2.y = quaternions[j*quaternion_offset+1];
		q2.z = quaternions[j*quaternion_offset+2];
		q2.w = quaternions[j*quaternion_offset+3];
		// calc the dot prod
		float dot_prod = dot(q1, q2);
		float4 pos_zero, neg;
		// avoiding branching, calc both
		pos_zero = q1 + weights[weight_index + j]*q2;
		neg = q1 - weights[weight_index + j]*q2;
		// simple ternary for assignment
		q1 = (dot_prod >= 0.0f) ? pos_zero : neg;
	}

	// normalize it
	q1 = normalize(q1);

	// convert it back to a rotation matrix
	/*

	this is row major, so we need column major
	Matrix<float, 4>(
    1.0f - 2.0f*qy*qy - 2.0f*qz*qz, 2.0f*qx*qy - 2.0f*qz*qw, 2.0f*qx*qz + 2.0f*qy*qw, 0.0f,
    2.0f*qx*qy + 2.0f*qz*qw, 1.0f - 2.0f*qx*qx - 2.0f*qz*qz, 2.0f*qy*qz - 2.0f*qx*qw, 0.0f,
    2.0f*qx*qz - 2.0f*qy*qw, 2.0f*qy*qz + 2.0f*qx*qw, 1.0f - 2.0f*qx*qx - 2.0f*qy*qy, 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f);

	column major
	Matrix<float, 4>(
    1.0f - 2.0f*qy*qy - 2.0f*qz*qz, 2.0f*qx*qy + 2.0f*qz*qw, 2.0f*qx*qz - 2.0f*qy*qw, 0.0f,
    2.0f*qx*qy - 2.0f*qz*qw, 1.0f - 2.0f*qx*qx - 2.0f*qz*qz, 2.0f*qy*qz + 2.0f*qx*qw, 0.0f,
    2.0f*qx*qz + 2.0f*qy*qw, 2.0f*qy*qz - 2.0f*qx*qw, 1.0f - 2.0f*qx*qx - 2.0f*qy*qy, 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f);
	
	*/

	float4 R[4];
	R[0] = (float4) (1.0f - 2.0f*q1.y*q1.y - 2.0f*q1.z*q1.z, \
		2.0f*q1.x*q1.y - 2.0f*q1.z*q1.w, \
		2.0f*q1.x*q1.z - 2.0f*q1.y*q1.w, \
		0.0f);
	R[1] = (float4) (2.0f*q1.x*q1.y - 2.0f*q1.z*q1.w,  \
		1.0f - 2.0f*q1.x*q1.x - 2.0f*q1.z*q1.z, \
		2.0f*q1.y*q1.z + 2.0f*q1.x*q1.w, \
		0.0f);
	R[2] = (float4) (2.0f*q1.x*q1.z + 2.0f*q1.y*q1.w,  \
		2.0f*q1.y*q1.z - 2.0f*q1.x*q1.w, \
		1.0f - 2.0f*q1.x*q1.x - 2.0f*q1.y*q1.y, \
		0.0f);
	R[3] = (float4) (0.0f, 0.0f, 0.0f, 1.0f);

	// perform the the weighted LBS sum to get R_prime and t_prime
	// setup the accumulator
	float4 R_prime_T_prime[4];
	for (unsigned int i = 0; i < 4; i++){
		R_prime_T_prime[i] = (float4) (0.0f);
	}
	// perform the weighted sum of the matricies
	for (unsigned int j = 0; j < numTransforms; j++){
		for (unsigned int i = 0; i < 4;  i++){
			R_prime_t_prime[i] += weights[weight_index+j]*matricies[j*4+i];
		}
	}

	// compute the translation (t) R_prime*CoR + t_prime - R*Cor
	// get the rotation only matrix from R_p_T_p
	float4 R_prime[4];
	for (unsigned int i = 0; i < 4; i++){
		R_prime[i] = R_prime_T_prime[i];
	}
	for (unsigned int i = 0; i < 3; i++){
		R_prime[i].w = 0.0f;
	}
	// get the translation from R_p_T_p
	float4 T_prime = (float4) (0.0f);
	T_prime.x = R_prime_T_prime[0].w;
	T_prime.y = R_prime_T_prime[1].w;
	T_prime.z = R_prime_T_prime[2].w;
	// compute R_prime * CoR
	float4 RpXCor = (float4) (0.0f);
	RpXCor.x = dot(CoRs[positionId], R_prime[0]);
	RpXCor.y = dot(CoRs[positionId], R_prime[1]);
	RpXCor.z = dot(CoRs[positionId], R_prime[2]);
	RpXCor.w = dot(CoRs[positionId], R_prime[3]);
	// compute R * CoR
	float4 RXCor = (float4) (0.0f);
	RXCor.x = dot(CoRs[positionId], R[0]);
	RXCor.y = dot(CoRs[positionId], R[1]);
	RXCor.z = dot(CoRs[positionId], R[2]);
	RXCor.w = dot(CoRs[positionId], R[3]);
	// compute t
	float4 T = (float4) (0.0f);
	T = RpXCor + T_prime - RXCor;
	// compute final pos, R*cur_pos + t
	finalPosition.x = dot(initialPosition, R[0]) + T.x;
	finalPosition.y = dot(initialPosition, R[1]) + T.y;
	finalPosition.z = dot(initialPosition, R[2]) + T.z;
	// put the final value into the output buffer
	finalPos[positionOffset] = finalPosition.x;
	finalPos[positionOffset+1] = finalPosition.y;
	finalPos[positionOffset+2] = finalPosition.z;

	// example code from offset project
	/*

	__global const float4* matrixInverse = &(matrices[4]);
	__global const float4* matrix = matrices;

	// point *= matrix inverse
	finalPosition.x = dot(initialPosition, matrixInverse[0]);
	finalPosition.y = dot(initialPosition, matrixInverse[1]);
	finalPosition.z = dot(initialPosition, matrixInverse[2]);

	// pt.y += weight
	finalPosition.y += weights[positionId];

	// point *= matrix
	// can't write back into finalPosition here, we need to use the same value to calculate xyz
	// instead write into global memory
	finalPos[positionOffset] = dot(finalPosition, matrix[0]);
	finalPos[positionOffset+1] = dot(finalPosition, matrix[1]);
	finalPos[positionOffset+2] = dot(finalPosition, matrix[2]);

	*/
}