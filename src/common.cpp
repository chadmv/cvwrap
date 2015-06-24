#include "common.h"

#include <maya/MGlobal.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnMesh.h>
#include <maya/MItMeshVertex.h>
#include <maya/MSelectionList.h>
#include <algorithm>
#include <cassert>
#include <set>
#include <queue>
#include <utility>

//#define NO_INTRINSICS
#define NORMALIZATION_INDEX -1

void StartProgress(const MString& title, unsigned int count) {
  if (MGlobal::mayaState() == MGlobal::kInteractive) {
    MString message = "progressBar -e -bp -ii true -st \"";
    message += title;
    message += "\" -max ";
    message += count;
    message += " $gMainProgressBar;";
    MGlobal::executeCommand(message);
  }
}


void StepProgress(int step) {
  if (MGlobal::mayaState() == MGlobal::kInteractive) {
    MString message = "progressBar -e -s ";
    message += step;
    message += " $gMainProgressBar;";
    MGlobal::executeCommand(message);
  }
}


bool ProgressCancelled() {
  if (MGlobal::mayaState() == MGlobal::kInteractive) {
    int cmdResult = 0;
    MGlobal::executeCommand("progressBar -query -isCancelled $gMainProgressBar", cmdResult);
    return cmdResult != 0;
  }
  return false;
}


void EndProgress() {
  if (MGlobal::mayaState() == MGlobal::kInteractive) {
    MGlobal::executeCommand("progressBar -e -ep $gMainProgressBar;");
  }
}


bool IsShapeNode(MDagPath& path) {
  return path.node().hasFn(MFn::kMesh) ||
         path.node().hasFn(MFn::kNurbsCurve) ||
         path.node().hasFn(MFn::kNurbsSurface);
}


MStatus GetShapeNode(MDagPath& path, bool intermediate) {
  MStatus status;

  if (IsShapeNode(path)) {
    // Start at the transform so we can honor the intermediate flag.
    path.pop();
  }

  if (path.hasFn(MFn::kTransform)) {
    unsigned int shapeCount;
    status = path.numberOfShapesDirectlyBelow(shapeCount);
    CHECK_MSTATUS_AND_RETURN_IT(status);

    for (unsigned int i = 0; i < shapeCount; ++i) {
      status = path.extendToShapeDirectlyBelow(i);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      // Make sure it is not an intermediate object.
      MFnDagNode fnNode(path, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      if ((!fnNode.isIntermediateObject() && !intermediate) || 
          (fnNode.isIntermediateObject() && intermediate)) {
        return MS::kSuccess;
      }
      // Go to the next shape.
      path.pop();
    }
  }

  // No valid shape node found.
  return MS::kFailure;
}


MStatus GetDagPath(MString& name, MDagPath& path) {
  MStatus status;
  MSelectionList list;
  status = MGlobal::getSelectionListByName(name, list);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = list.getDagPath(0, path);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  return MS::kSuccess;
}


MStatus DeleteIntermediateObjects(MDagPath& path) {
  MStatus status;
  if (IsShapeNode(path)) {
    path.pop();
  }
  if (path.hasFn(MFn::kTransform)) {
    unsigned int shapeCount;
    status = path.numberOfShapesDirectlyBelow(shapeCount);
    CHECK_MSTATUS_AND_RETURN_IT(status);

    for (unsigned int i = 0; i < shapeCount; ++i) {
      status = path.extendToShapeDirectlyBelow(i);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      MFnDagNode fnNode(path, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      path.pop();
      if (fnNode.isIntermediateObject()) {
        status = MGlobal::executeCommand("delete " + fnNode.partialPathName());
        CHECK_MSTATUS_AND_RETURN_IT(status);
      }
    }
  }
  return MS::kSuccess;
}


void GetBarycentricCoordinates(const MPoint& P, const MPoint& A, const MPoint& B, const MPoint& C,
                               BaryCoords& coords) {
  // Compute the normal of the triangle
  MVector N = (B - A) ^ (C - A);
  MVector unitN = N.normal();

  // Compute twice area of triangle ABC
  double areaABC = unitN * N;

  if (areaABC == 0.0) {
    // If the triangle is degenerate, just use one of the points.
    coords[0] = 1.0f;
    coords[1] = 0.0f;
    coords[2] = 0.0f;
    return;
  }

  // Compute a
  double areaPBC = unitN * ((B - P) ^ (C - P));
  coords[0] = (float)(areaPBC / areaABC);

  // Compute b
  double areaPCA = unitN * ((C - P) ^ (A - P));
  coords[1] = (float)(areaPCA / areaABC);

  // Compute c
  coords[2] = 1.0f - coords[0] - coords[1];
}


MStatus GetAdjacency(MDagPath& pathMesh, std::vector<std::set<int> >& adjacency) {
  MStatus status;
  // Get mesh adjacency.  The adjacency will be all vertex ids on the connected faces.
  MItMeshVertex itVert(pathMesh, MObject::kNullObj, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MFnMesh fnMesh(pathMesh, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  adjacency.resize(itVert.count());
  for (; !itVert.isDone(); itVert.next()) {
    MIntArray faces;
    status = itVert.getConnectedFaces(faces);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    adjacency[itVert.index()].clear();
    // Put the vertex ids in a set to avoid duplicates
    for (unsigned int j = 0; j < faces.length(); ++j) {
      MIntArray vertices;
      fnMesh.getPolygonVertices(faces[j], vertices);
      for (unsigned int k = 0; k < vertices.length(); ++k) {
        if (vertices[k] != itVert.index()) {
          adjacency[itVert.index()].insert(vertices[k]);
        }
      }
    }
  }
  return MS::kSuccess;
}


/**
  Used in the CrawlSurface function to keep track of where we are crawling.
*/
struct CrawlData {
  MPoint sourcePosition;  /**< Where the crawl iteration came from. */
  double crawlDistance;  /**< How far this crawl iteration has traveled. */
  int nextVertex;  /**< Where this crawl iteration should go next. */
};


MStatus CrawlSurface(const MPoint& startPoint, const MIntArray& vertexIndices, MPointArray& points, double maxDistance,
                     std::vector<std::set<int> >& adjacency, std::map<int, double>& distances) {
  MStatus status;
  distances[NORMALIZATION_INDEX] = 0.0; // -1 will represent our hit point.
  double minStartDistance = 999999.0;
  unsigned int minStartIndex = 0;

  // Instead of a recursive function, which can get pretty slow, we'll use a queue to keep
  // track of where we are going and where we are coming from.
  std::queue<CrawlData> verticesToVisit;
  // Add the initial crawl paths to the queue.
  for (unsigned int i = 0; i < vertexIndices.length(); ++i) {
    double distance = startPoint.distanceTo(points[vertexIndices[i]]);
    // Only crawl to the starting vertices if they are within the radius.
    if (distance <= maxDistance) {
      CrawlData root = {startPoint, distance, vertexIndices[i]};
      verticesToVisit.push(root);
      distances[vertexIndices[i]] = distance;
    }
    // Track the minimum start distance in case we need to add the closest vertex below.
    // The minimum must be greater than 0 to make sure we do not use the vertex that is the
    // same as the startPoint which would create an invalid up vector.
    if (distance < minStartDistance && distance > 0.000001) {
      minStartDistance = distance;
      minStartIndex = vertexIndices[i];
    }
  }
  // If we didn't even reach a vertex in the hit face, or the startPoint is equal to a vertex
  // on the face, add the closest vertex so we can calculate a proper up vector
  if (verticesToVisit.size() <= 1) {
    CrawlData root = {startPoint, maxDistance - 0.001, minStartIndex};
    verticesToVisit.push(root);
    distances[minStartIndex] = maxDistance - 0.001;
  }
  while (verticesToVisit.size()) {
    CrawlData next = verticesToVisit.front();
    verticesToVisit.pop();

    // Extract the data out of the crawl struct
    int idx = next.nextVertex;
    MPoint& pt = points[idx];
    MPoint sourcePoint = next.sourcePosition;
    double currentCrawlDistance = next.crawlDistance;

    currentCrawlDistance += sourcePoint.distanceTo(pt);
    if (currentCrawlDistance >= maxDistance) {
      // If this vertex is outside the radius, no need to crawl anymore from that vertex.
      continue;
    }
    double& savedDistance = distances[idx];
    if (currentCrawlDistance <= savedDistance || savedDistance == 0.0) {
      // If this current crawl distance is less then the distance we have saved for this
      // vertex, use this new crawl distance instead.
      savedDistance = currentCrawlDistance;
    } else {
      // A smaller distance is already stored so we don't want to crawl
      // from this vertex any further.
      continue;
    }
    // Crawl the adjacent vertices
    std::set<int>::iterator iter;
    for (iter = adjacency[idx].begin(); iter != adjacency[idx].end(); ++iter) {
      CrawlData data = {pt, currentCrawlDistance, *iter};
      verticesToVisit.push(data);
    }
  }
  assert(distances.size() > 0);
  
  return MS::kSuccess;
}

bool SampleSort(std::pair<int, double> lhs, std::pair<int, double> rhs) {
  // Ensure that the normalization sample comes last.
  return (lhs.second < rhs.second) || rhs.first == NORMALIZATION_INDEX; 
}

void CalculateSampleWeights(const std::map<int, double>& distances, double radius,
                            MIntArray& vertexIds, MDoubleArray& weights) {
  
  std::map<int, double>::const_iterator itDistance;
  std::vector<std::pair<int, double> > samples;
  for (itDistance = distances.begin();
        itDistance != distances.end();
        itDistance++) {
    double x = itDistance->second;
    double w = 1.0 - (x/radius);
    samples.push_back(std::pair<int, double>(itDistance->first, w));
  }

  // Make the samples a multiple of 4 so we can use fast intrinsics!
  int remainder = 4 - ((samples.size()-1) % 4);
  for (int i = 0; i < remainder; ++i) {
    samples.push_back(std::pair<int, double>(0, 0.0));
  }

  unsigned int length = (unsigned int)samples.size();
  weights.setLength(length);
  vertexIds.setLength(length);
  std::sort(samples.begin(), samples.end(), SampleSort);
  std::vector<std::pair<int, double> >::iterator iter;
  int ii = 0;
  double sum = 0.0;
  for (iter = samples.begin(); iter != samples.end(); ++iter, ++ii) {
    vertexIds[ii] = (*iter).first;
    weights[ii] = (*iter).second;
    sum += (*iter).second;
  }
  assert(sum > 0.0);
  // Normalize the weights
  for (unsigned int i = 0; i < weights.length(); ++i) {
    weights[i] /= sum;
  }
}


void CreateMatrix(const MPoint& origin, const MVector& normal, const MVector& up,
                  MMatrix& matrix) {
  const MPoint& t = origin;
  const MVector& y = normal;
  MVector x = y ^ up;
  MVector z = y ^ x;
  matrix[0][0] = x.x; matrix[0][1] = x.y; matrix[0][2] = x.z; matrix[0][3] = 0.0;
  matrix[1][0] = y.x; matrix[1][1] = y.y; matrix[1][2] = y.z; matrix[1][3] = 0.0;
  matrix[2][0] = z.x; matrix[2][1] = z.y; matrix[2][2] = z.z; matrix[2][3] = 0.0;
  matrix[3][0] = t.x; matrix[3][1] = t.y; matrix[3][2] = t.z; matrix[3][3] = 1.0;
}


void CalculateBasisComponents(const MDoubleArray& weights, const BaryCoords& coords,
                              const MIntArray& triangleVertices, const MPointArray& points,
                              const MFloatVectorArray& normals, const MIntArray& sampleIds,
                              double* alignedStorage,
                              MPoint& origin, MVector& up, MVector& normal) {
  // Start with the recreated point and normal using the barycentric coordinates of the hit point.
  unsigned int hitIndex = weights.length()-1;
#ifdef NO_INTRINSICS
  MPoint hitPoint;
  MVector hitNormal;
  for (int i = 0; i < 3; ++i) {
    hitPoint += points[triangleVertices[i]] * coords[i];
    hitNormal += MVector(normals[triangleVertices[i]]) * coords[i];
  }
  // Create the barycentric point and normal.
  origin = hitPoint * weights[hitIndex];
  normal = hitNormal * weights[hitIndex];
  // Then use the weighted adjacent data.
  for (unsigned int j = 0; j < hitIndex; j++) {
    origin += MVector(points[sampleIds[j]]) * weights[j];
    normal += MVector(normals[sampleIds[j]]) * weights[j];
  }

  // Calculate the up vector
  up = (hitPoint - origin) * weights[hitIndex];
  for (unsigned int j = 0; j < hitIndex; j++) {
    up += (points[sampleIds[j]] - origin) * weights[j];
  }
#else
  __m256d hitPointV = Dot4<MPoint>(coords[0], coords[1], coords[2], 0.0,
                                points[triangleVertices[0]], points[triangleVertices[1]],
                                points[triangleVertices[2]], MPoint::origin);
  __m256d hitNormalV = Dot4<MVector>(coords[0], coords[1], coords[2], 0.0,
                                normals[triangleVertices[0]], normals[triangleVertices[1]],
                                normals[triangleVertices[2]], MVector::zero);
  __m256d hitWeightV = _mm256_set1_pd(weights[hitIndex]);
  // Create the barycentric point and normal.
  __m256d originV = _mm256_mul_pd(hitPointV, hitWeightV);
  __m256d normalV = _mm256_mul_pd(hitNormalV, hitWeightV);
  // Then use the weighted adjacent data.
  for (unsigned int j = 0; j < hitIndex; j += 4) {
    __m256d tempOrigin = Dot4<MPoint>(weights[j], weights[j+1], weights[j+2], weights[j+3],
                                      points[sampleIds[j]], points[sampleIds[j+1]],
                                      points[sampleIds[j+2]], points[sampleIds[j+3]]);
    originV = _mm256_add_pd(tempOrigin, originV);

    __m256d tempNormal = Dot4<MVector>(weights[j], weights[j+1], weights[j+2], weights[j+3],
                                       normals[sampleIds[j]], normals[sampleIds[j+1]],
                                       normals[sampleIds[j+2]], normals[sampleIds[j+3]]);
    normalV = _mm256_add_pd(tempNormal, normalV);
  }

  //double values[4];
  _mm256_store_pd(alignedStorage, originV);
  origin.x = alignedStorage[0];
  origin.y = alignedStorage[1];
  origin.z = alignedStorage[2];
  _mm256_store_pd(alignedStorage, normalV);
  normal.x = alignedStorage[0];
  normal.y = alignedStorage[1];
  normal.z = alignedStorage[2];

  // Calculate the up vector
  __m256d upV = _mm256_mul_pd(_mm256_sub_pd(hitPointV, originV), hitWeightV);
  MVector v1, v2, v3, v4;
  for (unsigned int j = 0; j < hitIndex; j += 4) {
    v1 = points[sampleIds[j]] - origin;
    v2 = points[sampleIds[j+1]] - origin;
    v3 = points[sampleIds[j+2]] - origin;
    v4 = points[sampleIds[j+3]] - origin;
    __m256d tempUp = Dot4<MVector>(weights[j], weights[j+1], weights[j+2], weights[j+3],
                                   v1, v2, v3, v4);
    upV = _mm256_add_pd(tempUp, upV);
  }
  _mm256_store_pd(alignedStorage, upV);
  up.x = alignedStorage[0];
  up.y = alignedStorage[1];
  up.z = alignedStorage[2];
#endif

  GetValidUp(weights, points, sampleIds, origin, normal, up);
}


void GetValidUp(const MDoubleArray& weights, const MPointArray& points,
                const MIntArray& sampleIds, const MPoint& origin, const MVector& normal,
                MVector& up) {
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
}