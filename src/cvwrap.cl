/*
	offset kernels
*/

__kernel void cvwrap(
	__global float* finalPos,            //float3
	__global const float* initialPos,    //float3
	__global const float* driverPoints,  //float3
	__global const float* driverNormals, //float3
	__global const float* paintWeights,
	__global const int* sampleCounts,
	__global const int* sampleOffsets,
	__global const int* sampleIds,
	__global const float* sampleWeights,
	__global const int* triangleVerts,  //int3
	__global const float* baryCoords,  //float3
	__global const float4* bindMatrices,
	const uint positionCount)
{
	unsigned int positionId = get_global_id(0);				// access finalPos and initialPos using this value
	if (positionId >= positionCount) {
    // We create an execute unit for more indices then we have data for, just exit early if this guy if one of the extras
    return;					
  }
	unsigned int positionOffset = positionId * 3;				// Base positions are float3 when they come in here!
	

  // Start with the recreated point and normal using the barycentric coordinates of the hit point.
  float baryA = baryCoords[positionOffset];
  float baryB = baryCoords[positionOffset+1];
  float baryC = baryCoords[positionOffset+2];
  int triVertA = triangleVerts[positionOffset];
  int triVertB = triangleVerts[positionOffset+1];
  int triVertC = triangleVerts[positionOffset+2];
  float hitPointX = driverPoints[triVertA] * baryA +
                    driverPoints[triVertB] * baryB +
                    driverPoints[triVertC] * baryC;
  float hitPointY = driverPoints[triVertA+1] * baryA +
                    driverPoints[triVertB+1] * baryB +
                    driverPoints[triVertC+1] * baryC;
  float hitPointZ = driverPoints[triVertA+2] * baryA +
                    driverPoints[triVertB+2] * baryB +
                    driverPoints[triVertC+2] * baryC;
	float hitNormalX = driverNormals[triVertA] * baryA +
                     driverNormals[triVertB] * baryB +
                     driverNormals[triVertC] * baryC;
  float hitNormalY = driverNormals[triVertA+1] * baryA +
                     driverNormals[triVertB+1] * baryB +
                     driverNormals[triVertC+1] * baryC;
  float hitNormalZ = driverNormals[triVertA+2] * baryA +
                     driverNormals[triVertB+2] * baryB +
                     driverNormals[triVertC+2] * baryC;

  // Create the barycentric point and normal.
 // int hitIndex = sampleOffsets[positionId] + sampleCounts[positionId] - 1;
 // float hitWeight = sampleWeights[hitIndex];
 // float originX = hitPointX * hitWeight;
 // float originY = hitPointY * hitWeight;
 // float originZ = hitPointZ * hitWeight;
 // float normalX = hitNormalX * hitWeight;
 // float normalY = hitNormalY * hitWeight;
 // float normalZ = hitNormalZ * hitWeight;

 // // Then use the weighted adjacent data.
 // for (uint j = sampleOffsets[positionId]; j < hitIndex; j++) {
 //   float sw = sampleWeights[j];
 //   originX += driverPoints[sampleIds[j]*3]   * sw;
 //   originY += driverPoints[sampleIds[j]*3+1] * sw;
 //   originZ += driverPoints[sampleIds[j]*3+2] * sw;

 //   normalX += driverNormals[sampleIds[j]*3]   * sw;
 //   normalY += driverNormals[sampleIds[j]*3+1] * sw;
 //   normalZ += driverNormals[sampleIds[j]*3+2] * sw;
 // }

 // // Calculate the up vector
 // float upX = (hitPointX - originX) * hitWeight;
 // float upY = (hitPointY - originY) * hitWeight;
 // float upZ = (hitPointZ - originZ) * hitWeight;
 // for (uint j = sampleOffsets[positionId]; j < hitIndex; j++) {
 //   float sw = sampleWeights[j];
 //   upX += (driverPoints[sampleIds[j]*3] - originX) * sw;
 //   upY += (driverPoints[sampleIds[j]*3+1] - originY) * sw;
 //   upZ += (driverPoints[sampleIds[j]*3+2] - originZ) * sw;
 // }

 // // Use float3 so we can use the built-in functions.  We are mostly using single floats
 // // because the preferred vector width of most gpu's these days is 1.
 // float3 up = (float3)(upX, upY, upZ);
 // float3 normal = (float3)(normalX, normalY, normalZ);
 // float3 unitUp = fast_normalize(up);
 // float upLength = fast_length(up);
 // if (fabs(dot(unitUp, normal) - 1.0f) > 0.001f && upLength > 0.0001f) {
 //   for (uint j = sampleOffsets[positionId]; j < hitIndex; j++) {
 //     up.x -= (driverPoints[sampleIds[j]*3] - originX) * sampleWeights[j];
 //     up.y -= (driverPoints[sampleIds[j]*3+1] - originY) * sampleWeights[j];
 //     up.z -= (driverPoints[sampleIds[j]*3+2] - originZ) * sampleWeights[j];
 //     unitUp = fast_normalize(up);
 //     upLength = fast_length(up);
 //     if (fabs(dot(unitUp, normal) - 1.0f) > 0.001f && upLength > 0.0001f) {
 //       // If the up and normal vectors are no longer parallel and the up vector has a length,
 //       // then we are good to go.
 //       break;
 //     }
 //   }
 //   up = fast_normalize(up);
 // } else {
 //   up = unitUp;
 // }

 // // Create the transform matrix
 // // Store by columns so we can use dot to multiply with the bind matrix
 // float3 x = cross(normal, up);
 // float3 z = cross(normal, x);
 // float4 matrix0 = (float4)(x.x, normal.x, z.x, originX);
 // float4 matrix1 = (float4)(x.y, normal.y, z.y, originY);
 // float4 matrix2 = (float4)(x.z, normal.z, z.z, originZ);
 // float4 matrix3 = (float4)(0.0f, 0.0f, 0.0f, 1.0f);

 // // TODO: scale matrix mult

 // // Multiply bindMatrix with matrix
	//__global const float4* bindMatrix = &(bindMatrices[positionId*4]);
 // float4 bindMatrix0 = bindMatrix[0];
 // float4 bindMatrix1 = bindMatrix[1];
 // float4 bindMatrix2 = bindMatrix[2];
 // float4 m0 = (float4)(dot(bindMatrix0, matrix0),
 //                      dot(bindMatrix0, matrix1),
 //                      dot(bindMatrix0, matrix2),
 //                      dot(bindMatrix0, matrix3));
 // float4 m1 = (float4)(dot(bindMatrix1, matrix0),
 //                      dot(bindMatrix1, matrix1),
 //                      dot(bindMatrix1, matrix2),
 //                      dot(bindMatrix1, matrix3));
 // float4 m2 = (float4)(dot(bindMatrix2, matrix0),
 //                      dot(bindMatrix2, matrix1),
 //                      dot(bindMatrix2, matrix2),
 //                      dot(bindMatrix2, matrix3));

	float4 initialPosition = (float4)(initialPos[positionOffset],
	                                  initialPos[positionOffset+1],
	                                  initialPos[positionOffset+2],
                                    1.0f);
	/*finalPos[positionOffset] = dot(initialPosition, m0);
	finalPos[positionOffset+1] = dot(initialPosition, m1);
	finalPos[positionOffset+2] = dot(initialPosition, m2);*/
  finalPos[positionOffset] = hitPointX;
	finalPos[positionOffset+1] = hitPointY;
	finalPos[positionOffset+2] = hitPointZ;
}