#include "common.h"

#include <maya/MGlobal.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnMesh.h>
#include <maya/MItMeshVertex.h>
#include <maya/MSelectionList.h>
#include <cassert>
#include <set>
#include <queue>
#include <utility>

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


double DistanceSquared(const MPoint& p1, const MPoint& p2) {
  double xx = p2.x - p1.x;
  double yy = p2.y - p1.y;
  double zz = p2.z - p1.z;
  return (xx*xx) + (yy*yy) + (zz*zz);
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
  // Instead of a recursive function, which can get pretty slow, we'll use a queue to keep
  // track of where we are going and where we are coming from.
  std::queue<CrawlData> verticesToVisit;
  // Add the initial crawl paths to the queue.
  for (unsigned int i = 0; i < vertexIndices.length(); ++i) {
    CrawlData root = {startPoint, 0.0, vertexIndices[i]};
    verticesToVisit.push(root);
    // Store the initial distances so we hit all the starting vertices.
    distances[vertexIndices[i]] = startPoint.distanceTo(points[vertexIndices[i]]);
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


void CalculateSampleWeights(const std::map<int, double>& distances, double radius,
                            MIntArray& vertexIds, MDoubleArray& weights) {
  unsigned int length = (unsigned int)distances.size();
  weights.setLength(length);
  vertexIds.setLength(length);
  int ii = 0;
  // Use a gaussian function to calculate the sample weights.
  // sigmaSquared corresponds to the width of the bell curve.
  double sigmaSquared = (radius * radius) * 0.75;
  std::map<int, double>::const_iterator itDistance;
  for (itDistance = distances.begin();
        itDistance != distances.end();
        itDistance++, ii++) {
    vertexIds[ii] = itDistance->first;
    double x = itDistance->second;
    weights[ii] = exp(-(x*x)/(sigmaSquared));
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
