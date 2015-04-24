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
#include <xmmintrin.h>
#include <emmintrin.h>
#include <immintrin.h>

#define NORMALIZATION_INDEX -1
#define MM256_SHUFFLE(fp3,fp2,fp1,fp0) (((fp3) << 3) | ((fp2) << 2) | \
                                     ((fp1) << 1) | ((fp0)))
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


MStatus GetShapeNode(MDagPath& path) {
  MStatus status;

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
      if (!fnNode.isIntermediateObject()) {
        return MS::kSuccess;
      }
      // Go to the next shape if this is an intermediate shape.
      path.pop();
    }
  } else if (IsShapeNode(path)) {
    // This is already a shape node.
    return MS::kSuccess;
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
    if (distance < minStartDistance) {
      minStartDistance = distance;
      minStartIndex = vertexIndices[i];
    }
  }
  // If we didn't even reach a vertex in the hit face, add the closest vertex so we can calculate
  // a proper up vector
  if (verticesToVisit.size() == 0) {
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

  unsigned int length = (unsigned int)distances.size();
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

__m256d Dot(const __m256d v1, const __m256d v2) {
  // dot = (w1*w2, z1*z2, y1*y2, x1*x2)
  __m256d dot = _mm256_mul_pd(v1, v2);
  // temp = (y1*y2, y1*y2, z1*z2, z1*z2)
  __m256d temp = _mm256_shuffle_pd(dot, dot, MM256_SHUFFLE(0, 0, 1, 1));
  // dot = (w1*w2, z1*z2, y1*y2, x1*x2+y1*y2)
  dot = _mm256_add_pd(dot, temp);
  // temp = (z1*z2, z1*z2, z1*z2, z1*z2)
  temp = _mm256_shuffle_pd(temp, temp, MM256_SHUFFLE(1, 1, 1, 1));
  // dot = (w1*w2, z1*z2, y1*y2, x1*x2+y1*y2+z1*z2)
  dot = _mm256_add_pd(dot, temp);
  __m256d result;
  // values = (x1*x2+y1*y2+z1*z2, x1*x2+y1*y2+z1*z2, x1*x2+y1*y2+z1*z2, x1*x2+y1*y2+z1*z2)
  result = _mm256_shuffle_pd(dot, dot, MM256_SHUFFLE(0, 0, 0, 0));
  return result;
}

template <typename T>
void CalculateBaryPointAVX(const MIntArray& triangleVertices, const BaryCoords& coords,
                           const T& p1, const T& p2, const T& p3,
                           double& x, double& y, double& z) {
  __m256d xxx = _mm256_set_pd(p1.x, p2.x, p3.x, 0.0);
  __m256d yyy = _mm256_set_pd(p1.y, p2.y, p3.y, 0.0);
  __m256d zzz = _mm256_set_pd(p1.z, p2.z, p3.z, 0.0);
  __m256d www = _mm256_set_pd(coords[0], coords[1], coords[2], 0.0);
  __m256d xy0 = _mm256_mul_pd(xxx, www);
  __m256d xy1 = _mm256_mul_pd(yyy, www);
  __m256d xy2 = _mm256_mul_pd(zzz, www);
  __m256d xy3 = _mm256_mul_pd(www, www); // Dummy
  // low to high: xy00+xy01 xy10+xy11 xy02+xy03 xy12+xy13
  __m256d temp01 = _mm256_hadd_pd(xy0, xy1);   
  // low to high: xy20+xy21 xy30+xy31 xy22+xy23 xy32+xy33
  __m256d temp23 = _mm256_hadd_pd(xy2, xy3);
  // low to high: xy02+xy03 xy12+xy13 xy20+xy21 xy30+xy31
  __m256d swapped = _mm256_permute2f128_pd(temp01, temp23, 0x21);
  // low to high: xy00+xy01 xy10+xy11 xy22+xy23 xy32+xy33
  __m256d blended = _mm256_blend_pd(temp01, temp23, 0x0C);
  __m256d dotproduct = _mm256_add_pd(swapped, blended);
  double values[4];
  _mm256_store_pd(values, dotproduct);
  x = values[0];
  y = values[1];
  z = values[2];
}

void CalculateBasisComponents(const MDoubleArray& weights, const BaryCoords& coords,
                              const MIntArray& triangleVertices, const MPointArray& points,
                              const MFloatVectorArray& normals, const MIntArray& sampleIds,
                              MPoint& origin, MVector& up, MVector& normal) {
  // Start with the recreated point and normal using the barycentric coordinates of the hit point.
  MPoint hitPoint;
  MVector hitNormal;
  CalculateBaryPointAVX<MPoint>(triangleVertices, coords,
                                points[triangleVertices[0]], points[triangleVertices[1]],
                                points[triangleVertices[2]], hitPoint.x, hitPoint.y, hitPoint.z);
  CalculateBaryPointAVX<MVector>(triangleVertices, coords,
                                normals[triangleVertices[0]], normals[triangleVertices[1]],
                                normals[triangleVertices[2]], hitNormal.x, hitNormal.y, hitNormal.z);
 

  for (int i = 0; i < 3; ++i) {
    //hitPoint += points[triangleVertices[i]] * coords[i];
    //hitNormal += MVector(normals[triangleVertices[i]]) * coords[i];
  }
  // Then use the weighted adjacent data.
  unsigned int hitIndex = weights.length()-1;
  for (unsigned int j = 0; j < hitIndex; j++) {
    normal += MVector(normals[sampleIds[j]]) * weights[j];
    origin += MVector(points[sampleIds[j]]) * weights[j];
  }
  // Add the recreated barycentric point and normal.
  normal += hitNormal * weights[hitIndex];
  origin += hitPoint * weights[hitIndex];

  for (unsigned int j = 0; j < hitIndex; j++) {
    up += (points[sampleIds[j]] - origin) * weights[j];
  }
  up += (hitPoint - origin) * weights[hitIndex];
  normal.normalize();

  GetValidUpAndNormal(weights, points, sampleIds, origin, up, normal);
}


void GetValidUpAndNormal(const MDoubleArray& weights, const MPointArray& points,
                         const MIntArray& sampleIds, const MPoint& origin, MVector& up,
                         MVector& normal) {
  MVector unitUp = up.normal();
  // Adjust up if it's parallel to normal or if it's zero length
  if (unitUp * normal == 1.0 || up.length() < 0.0001) {
    for (unsigned int j = 0; j < weights.length()-1; ++j) {
      if (unitUp * normal != 1.0 && up.length() > 0.0001) {
        // If the up and normal vectors are no longer parallel and the up vector has a length,
        // then we are good to go.
        break;
      }
      up -= (points[sampleIds[j]] - origin) * weights[j];
      unitUp = up.normal();
    }
    up.normalize();
  } else {
    up = unitUp;
  }
}