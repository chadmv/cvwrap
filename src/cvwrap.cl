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
    MPoint hitPoint;
    MVector hitNormal;
    for (int i = 0; i < 3; ++i) {
      hitPoint += points[triangleVertices[i]] * coords[i];
      hitNormal += MVector(normals[triangleVertices[i]]) * coords[i];
    }
  */
  float baryA = baryCoords[positionOffset];
  float baryB = baryCoords[positionOffset+1];
  float baryC = baryCoords[positionOffset+2];
  int triVertA = triangleVerts[positionOffset] * 3;
  int triVertB = triangleVerts[positionOffset+1] * 3;
  int triVertC = triangleVerts[positionOffset+2] * 3;
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
  /*
    Equivalent CPU code:
    ====================
    unsigned int hitIndex = weights.length()-1;
    origin = hitPoint * weights[hitIndex];
    normal = hitNormal * weights[hitIndex];
  */
  int offset = sampleOffsets[positionId];
  int hitIndex = offset + sampleCounts[positionId] - 1;
  float hitWeight = sampleWeights[hitIndex];
  float originX = hitPointX * hitWeight;
  float originY = hitPointY * hitWeight;
  float originZ = hitPointZ * hitWeight;
  float normalX = hitNormalX * hitWeight;
  float normalY = hitNormalY * hitWeight;
  float normalZ = hitNormalZ * hitWeight;

  // Then use the weighted adjacent data.
  /*
    Equivalent CPU code:
    ====================
    for (unsigned int j = 0; j < hitIndex; j++) {
      origin += MVector(points[sampleIds[j]]) * weights[j];
      normal += MVector(normals[sampleIds[j]]) * weights[j];
    }
  */
  for (int j = offset; j < hitIndex; j++) {
    float sw = sampleWeights[j];
    int sampleId = sampleIds[j] * 3;
    originX += driverPoints[sampleId]   * sw;
    originY += driverPoints[sampleId+1] * sw;
    originZ += driverPoints[sampleId+2] * sw;

    normalX += driverNormals[sampleId]   * sw;
    normalY += driverNormals[sampleId+1] * sw;
    normalZ += driverNormals[sampleId+2] * sw;
  }

  // Calculate the up vector
  /*
    Equivalent CPU code:
    ====================
    up = (hitPoint - origin) * weights[hitIndex];
    for (unsigned int j = 0; j < hitIndex; j++) {
      up += (points[sampleIds[j]] - origin) * weights[j];
    }
  */
  float upX = (hitPointX - originX) * hitWeight;
  float upY = (hitPointY - originY) * hitWeight;
  float upZ = (hitPointZ - originZ) * hitWeight;
  for (int j = offset; j < hitIndex; j++) {
    float sw = sampleWeights[j];
    int sampleId = sampleIds[j] * 3;
    upX += (driverPoints[sampleId] - originX) * sw;
    upY += (driverPoints[sampleId+1] - originY) * sw;
    upZ += (driverPoints[sampleId+2] - originZ) * sw;
  }

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
  // Store by columns so we can use dot to multiply with the bind matrix
  float3 x = cross(normal, up);
  float3 z = cross(normal, x);
  float4 matrix0 = (float4)(x.x, normal.x, z.x, originX);
  float4 matrix1 = (float4)(x.y, normal.y, z.y, originY);
  float4 matrix2 = (float4)(x.z, normal.z, z.z, originZ);
  float4 matrix3 = (float4)(0.0f, 0.0f, 0.0f, 1.0f);

 // // TODO: scale matrix mult

  // Multiply bindMatrix with matrix
  /*
    Equivalent CPU code:
    ====================
    MPoint newPt = (points[i] * (bindMatrices[index] * matrix));
  */
  float4 bindMatrix0 = bindMatrices[positionId*4];
  float4 bindMatrix1 = bindMatrices[positionId*4+1];
  float4 bindMatrix2 = bindMatrices[positionId*4+2];
  float4 bindMatrix3 = bindMatrices[positionId*4+3];
  float4 m0 = (float4)(dot(bindMatrix0, matrix0),
                       dot(bindMatrix0, matrix1),
                       dot(bindMatrix0, matrix2),
                       dot(bindMatrix0, matrix3));
  float4 m1 = (float4)(dot(bindMatrix1, matrix0),
                       dot(bindMatrix1, matrix1),
                       dot(bindMatrix1, matrix2),
                       dot(bindMatrix1, matrix3));
  float4 m2 = (float4)(dot(bindMatrix2, matrix0),
                       dot(bindMatrix2, matrix1),
                       dot(bindMatrix2, matrix2),
                       dot(bindMatrix2, matrix3));
  float4 m3 = (float4)(dot(bindMatrix3, matrix0),
                       dot(bindMatrix3, matrix1),
                       dot(bindMatrix3, matrix2),
                       dot(bindMatrix3, matrix3));

  float4 initialPosition = (float4)(initialPos[positionOffset],
                                    initialPos[positionOffset+1],
                                    initialPos[positionOffset+2],
                                    1.0f);
  finalPos[positionOffset] = dot(initialPosition, (float4)(m0.x, m1.x, m2.x, m3.x));
  finalPos[positionOffset+1] = dot(initialPosition, (float4)(m0.y, m1.y, m2.y, m3.y));
  finalPos[positionOffset+2] = dot(initialPosition, (float4)(m0.z, m1.z, m2.z, m3.z));
}