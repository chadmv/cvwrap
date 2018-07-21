/*
  cvwrap kernel
*/

__kernel void cvwrap(__global float* finalPos,
                     __global const float* initialPos,
                     __global const float* driverPoints,
                     __global const float* driverNormals,
                     __global const float* paintWeights,
                     __global const int* sampleCounts,
                     __global const int* sampleOffsets,
                     __global const int* sampleIds,
                     __global const float* sampleWeights,
                     __global const int* triangleVerts,
                     __global const float* baryCoords,
                     __global const float4* bindMatrices,
                     __global const float4* drivenMatrices,
                     const float envelope,
                     const uint positionCount) {
  unsigned int positionId = get_global_id(0);
  if (positionId >= positionCount) {
    return;          
  }
  unsigned int positionOffset = positionId * 3;

  // Start with the recreated point and normal using the barycentric coordinates of the hit point.
  /*
    Equivalent CPU code:
    ====================
    MVector hitNormal;
    for (int i = 0; i < 3; ++i) {
      origin += points[triangleVertices[i]] * coords[i];
      hitNormal += MVector(normals[triangleVertices[i]]) * coords[i];
    }
  */
  float baryA = baryCoords[positionOffset];
  float baryB = baryCoords[positionOffset+1];
  float baryC = baryCoords[positionOffset+2];
  int triVertA = triangleVerts[positionOffset] * 3;
  int triVertB = triangleVerts[positionOffset+1] * 3;
  int triVertC = triangleVerts[positionOffset+2] * 3;
  float originX = driverPoints[triVertA] * baryA +
                  driverPoints[triVertB] * baryB +
                  driverPoints[triVertC] * baryC;
  float originY = driverPoints[triVertA+1] * baryA +
                  driverPoints[triVertB+1] * baryB +
                  driverPoints[triVertC+1] * baryC;
  float originZ = driverPoints[triVertA+2] * baryA +
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

  /*
    Equivalent CPU code:
    ====================
    unsigned int hitIndex = weights.length()-1;
    normal = hitNormal * weights[hitIndex];
  */
  int offset = sampleOffsets[positionId];
  int hitIndex = offset + sampleCounts[positionId] - 1;
  float hitWeight = sampleWeights[hitIndex];
  float normalX = hitNormalX * hitWeight;
  float normalY = hitNormalY * hitWeight;
  float normalZ = hitNormalZ * hitWeight;

  // Use crawl data to calculate normal
  /*
    Equivalent CPU code:
    ====================
    for (unsigned int j = 0; j < hitIndex; j++) {
      normal += MVector(normals[sampleIds[j]]) * weights[j];
    }
  */
  for (int j = offset; j < hitIndex; j++) {
    float sw = sampleWeights[j];
    int sampleId = sampleIds[j] * 3;
    normalX += driverNormals[sampleId]   * sw;
    normalY += driverNormals[sampleId+1] * sw;
    normalZ += driverNormals[sampleId+2] * sw;
  }

  // Calculate the up vector
  /*
    Equivalent CPU code:
    ====================
    up = ((points[triangleVertices[0]] + points[triangleVertices[1]]) * 0.5) - origin;
  */
  float upX = ((driverPoints[triVertA] + driverPoints[triVertB]) * 0.5f) - originX;
  float upY = ((driverPoints[triVertA+1] + driverPoints[triVertB+1]) * 0.5f) - originY;
  float upZ = ((driverPoints[triVertA+2] + driverPoints[triVertB+2]) * 0.5f) - originZ;

  // Use float3 so we can use the built-in functions.  We are mostly using single floats
  // because the preferred vector width of most gpu's these days is 1.
  /*
    Equivalent CPU code:
    ====================
    MVector unitUp = up.normal();
    // Adjust up if it's parallel to normal or if it's zero length
    if (abs((unitUp * normal) - 1.0) < 0.001 || up.length() < 0.0001) {
      for (unsigned int j = 0; j < weights.length()-1; ++j) {
        up -= (points[sampleIds[j]] - origin) * weights[j];
        unitUp = up.normal();
        if (abs((unitUp * normal) - 1.0) > 0.001 && up.length() > 0.0001) {
          // If the up and normal vectors are no longer parallel and the up vector has a length,
          // then we are good to go.
          break;
        }
      }
      up.normalize();
    } else {
      up = unitUp;
    }
  */
  float3 up = (float3)(upX, upY, upZ);
  float3 normal = (float3)(normalX, normalY, normalZ);
  normal = normalize(normal);
  float3 unitUp = normalize(up);
  float upLength = length(up);
  if (fabs(dot(unitUp, normal) - 1.0f) < 0.001f || upLength < 0.0001f) {
    for (int j = offset; j < hitIndex; j++) {
      float sw = sampleWeights[j];
      int sampleId = sampleIds[j] * 3;
      up.x -= (driverPoints[sampleId] - originX) * sw;
      up.y -= (driverPoints[sampleId+1] - originY) * sw;
      up.z -= (driverPoints[sampleId+2] - originZ) * sw;
      unitUp = normalize(up);
      upLength = length(up);
      if (fabs(dot(unitUp, normal) - 1.0f) > 0.001f && upLength > 0.0001f) {
        // If the up and normal vectors are no longer parallel and the up vector has a length,
        // then we are good to go.
        break;
      }
    }
    up = normalize(up);
  } else {
    up = unitUp;
  }

  // Create the transform matrix
  // Store by columns so we can use dot to multiply with the scale matrix
  float3 x = cross(normal, up);
  float3 z = cross(x, normal);
  x = normalize(x);
  z = normalize(z);

  float4 matrix0 = (float4)(x.x, normal.x, z.x, originX);
  float4 matrix1 = (float4)(x.y, normal.y, z.y, originY);
  float4 matrix2 = (float4)(x.z, normal.z, z.z, originZ);
  float4 matrix3 = (float4)(0.0f, 0.0f, 0.0f, 1.0f);

  // Scale matrix mult
  /*
    Equivalent CPU code:
    ====================
    matrix = scaleMatrix * matrix;
  */
  __global const float4* scaleMatrix = &(drivenMatrices[8]);
  float4 scaleMatrix0 = (float4)(dot(scaleMatrix[0], matrix0),
                        dot(scaleMatrix[0], matrix1),
                        dot(scaleMatrix[0], matrix2),
                        dot(scaleMatrix[0], matrix3));
  float4 scaleMatrix1 = (float4)(dot(scaleMatrix[1], matrix0),
                        dot(scaleMatrix[1], matrix1),
                        dot(scaleMatrix[1], matrix2),
                        dot(scaleMatrix[1], matrix3));
  float4 scaleMatrix2 = (float4)(dot(scaleMatrix[2], matrix0),
                        dot(scaleMatrix[2], matrix1),
                        dot(scaleMatrix[2], matrix2),
                        dot(scaleMatrix[2], matrix3));
  float4 scaleMatrix3 = (float4)(dot(scaleMatrix[3], matrix0),
                        dot(scaleMatrix[3], matrix1),
                        dot(scaleMatrix[3], matrix2),
                        dot(scaleMatrix[3], matrix3));
  // Transpose so we can dot with bindMatrices
  float4 smX = (float4)(scaleMatrix0.x, scaleMatrix1.x, scaleMatrix2.x, scaleMatrix3.x);
  float4 smY = (float4)(scaleMatrix0.y, scaleMatrix1.y, scaleMatrix2.y, scaleMatrix3.y);
  float4 smZ = (float4)(scaleMatrix0.z, scaleMatrix1.z, scaleMatrix2.z, scaleMatrix3.z);
  float4 smW = (float4)(scaleMatrix0.w, scaleMatrix1.w, scaleMatrix2.w, scaleMatrix3.w);

  // Multiply bindMatrix with matrix
  /*
    Equivalent CPU code:
    ====================
    MPoint newPt = ((points[i]  * drivenMatrix) * (bindMatrices[index] * matrix)) * drivenInverseMatrix;
  */
  float4 bm0 = bindMatrices[positionId*4];
  float4 bm1 = bindMatrices[positionId*4+1];
  float4 bm2 = bindMatrices[positionId*4+2];
  float4 bm3 = bindMatrices[positionId*4+3];
  float4 m0 = (float4)(dot(bm0, smX), dot(bm0, smY), dot(bm0, smZ), dot(bm0, smW));
  float4 m1 = (float4)(dot(bm1, smX), dot(bm1, smY), dot(bm1, smZ), dot(bm1, smW));
  float4 m2 = (float4)(dot(bm2, smX), dot(bm2, smY), dot(bm2, smZ), dot(bm2, smW));
  float4 m3 = (float4)(dot(bm3, smX), dot(bm3, smY), dot(bm3, smZ), dot(bm3, smW));

  float4 initialPosition = (float4)(initialPos[positionOffset],
                                    initialPos[positionOffset+1],
                                    initialPos[positionOffset+2],
                                    1.0f);
  __global const float4* drivenInverseMatrix = &(drivenMatrices[4]);
	__global const float4* drivenMatrix = drivenMatrices;
  float4 worldPt = (float4)(dot(initialPosition, drivenMatrix[0]),
                            dot(initialPosition, drivenMatrix[1]),
                            dot(initialPosition, drivenMatrix[2]),
                            dot(initialPosition, drivenMatrix[3]));
  worldPt = (float4)(dot(worldPt, (float4)(m0.x, m1.x, m2.x, m3.x)),
                     dot(worldPt, (float4)(m0.y, m1.y, m2.y, m3.y)),
                     dot(worldPt, (float4)(m0.z, m1.z, m2.z, m3.z)),
                     dot(worldPt, (float4)(m0.w, m1.w, m2.w, m3.w)));
  float3 newPt = (float3)(dot(worldPt, drivenInverseMatrix[0]),
                          dot(worldPt, drivenInverseMatrix[1]),
                          dot(worldPt, drivenInverseMatrix[2]));
  /*
    Equivalent CPU code:
    ====================
    points[i] = points[i] + ((newPt - points[i]) * paintWeights[i] * env);
  */
  float weight = paintWeights[positionId] * envelope;
  finalPos[positionOffset] = initialPosition.x + ((newPt.x - initialPosition.x) * weight);
  finalPos[positionOffset+1] = initialPosition.y + ((newPt.y - initialPosition.y) * weight);
  finalPos[positionOffset+2] = initialPosition.z + ((newPt.z - initialPosition.z) * weight);
}