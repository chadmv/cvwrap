/**
  Contains various helper functions.
*/

#ifndef CVWRAP_COMMON_H
#define CVWRAP_COMMON_H

#include <maya/MDagPath.h>
#include <maya/MDoubleArray.h>
#include <maya/MIntArray.h>
#include <maya/MMatrix.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MString.h>
#include <map>
#include <vector>
#include <set>
/**
  Helper function to start a new progress bar.
  @param[in] title Status title.
  @param[in] count Progress bar maximum count.
*/
void StartProgress(const MString& title, unsigned int count);


/**
  Helper function to increase the progress bar by the specified amount.
  @param[in] step Step amount.
*/
void StepProgress(int step);


/**
  Check if the progress has been cancelled.
  @return true if the progress has been cancelled.
*/
bool ProgressCancelled();


/**
  Ends any running progress bar.
*/
void EndProgress();


/**
  Checks if the path points to a shape node.
  @param[in] path A dag path.
  @return true if the path points to a shape node.
 */
bool IsShapeNode(MDagPath& path);


/**
  Ensures that the given dag path points to a non-intermediate shape node.
  @param[in,out] path Path to a dag node that could be a transform or a shape.
  On return, the path will be to a shape node if one exists.
  @return MStatus.
 */
MStatus GetShapeNode(MDagPath& path);


/**
  Get the MDagPath of an object.
  @param[in] name Name of a dag node.
  @param[out] path Storage for the dag path.
 */
MStatus GetDagPath(MString& name, MDagPath& path);


/**
  Delete all intermediate shapes of the given dag path.
  @param[in] path MDagPath.
 */
MStatus DeleteIntermediateObjects(MDagPath& path);

double DistanceSquared(const MPoint& p1, const MPoint& p2);

/**
  Helper struct to hold the 3 barycentric coordinates.
*/
struct BaryCoords {
  float coords[3];
  float operator[](int index) const { return coords[index]; }
  float& operator[](int index) { return coords[index]; }
};

/**
  Get the barycentric coordinates of point P in the triangle specified by points A,B,C.
  @param[in] P The sample point.
  @param[in] A Triangle point.
  @param[in] B Triangle point.
  @param[in] C Triangle point.
  @param[out] coords Barycentric coordinates storage.
*/
void GetBarycentricCoordinates(const MPoint& P, const MPoint& A, const MPoint& B, const MPoint& C,
                               BaryCoords& coords);

/**
  Get the vertex adjacency of the specified mesh.  The vertex adjacency are the vertex ids
  of the connected faces of each vertex.
  @param[in] pathMesh Path to a mesh.
  @param[out] adjacency Ajdancency storage of the adjancency per vertex id.
 */
MStatus GetAdjacency(MDagPath& pathMesh, std::vector<std::set<int> >& adjacency);


/**
  Crawls the surface to find all the points within the sampleradius.
  @param[in] startPoint The position from which to start the crawl.
  @param[in] vertexIndices The starting vertex indices we want to crawl to.
  @param[in] points The array of all the mesh points.
  @param[in] maxDistance The maximum crawl distance.
  @param[in] adjacency Vertex adjacency data from the GetAdjacency function.
  @param[out] distances Storage for the distances to the crawled points.
  @return MStatus
 */
MStatus CrawlSurface(const MPoint& startPoint, const MIntArray& vertexIndices, MPointArray& points, double maxDistance,
                     std::vector<std::set<int> >& adjacency, std::map<int, double>& distances);


/**
  Calculates a weight for each vertex within the crawl sample radius.  Vertices that are further
  away from the origin should have a lesser effect than vertices closer to the origin.
  @param[in] distances Crawl distances calculated from CrawlSurface.
  @param[in] radius Sample radius.
  @param[out] vertexIds Storage for the vertex ids sampled during the crawl.
  @param[out] weights Storage for the calculated weights of each sampled vertex.
*/
 void CalculateSampleWeights(const std::map<int, double>& distances, double radius,
                             MIntArray& vertexIds, MDoubleArray& weights);

 /**
   Creates an orthonormal basis using the given point and two axes.
   @param[in] origin Position.
   @param[in] normal Normal vector.
   @param[in] up Up vector.
   @param[out] matrix Generated matrix.
 */
 void CreateMatrix(const MPoint& origin, const MVector& normal, const MVector& up,
                   MMatrix& matrix);

#endif