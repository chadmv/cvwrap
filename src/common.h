/**
  Contains various helper functions.
*/

#ifndef CVWRAP_COMMON_H
#define CVWRAP_COMMON_H

#include <maya/MDagPath.h>
#include <maya/MDoubleArray.h>
#include <maya/MFloatVectorArray.h>
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

/**
  Calculates the components necessary to create a wrap basis matrix.
  @param[in] weights The sample weights array from the wrap binding.
  @param[in] points The driver point array.
  @param[in] normals The driver per-vertex normal array.
  @param[in] sampleIds The vertex ids on the driver of the current sample.
  @param[out] origin The origin of the coordinate system.
  @param[out] up The up vector of the coordinate system.
  @param[out] normal The normal vector of the coordinate system.
*/
void CalculateBasisComponents(const MDoubleArray& weights, const MPointArray& points,
                              const MFloatVectorArray& normals, const MIntArray& sampleIds,
                              MPoint& origin, MVector& up, MVector& normal);

/**
  Ensures that the up and normal vectors are perpendicular to each other.
  @param[in] weights The sample weights array from the wrap binding.
  @param[in] points The driver point array.
  @param[in] sampleIds The vertex ids on the driver of the current sample.
  @param[in] origin The origin of the coordinate system.
  @param[out] up The up vector of the coordinate system.
  @param[out] normal The normal vector of the coordinate system.
*/
void GetValidUpAndNormal(const MDoubleArray& weights, const MPointArray& points,
                         const MIntArray& sampleIds, const MPoint& origin, MVector& up,
                         MVector& normal);


template <typename T>
struct ThreadData {
  unsigned int start;
  unsigned int end;
  unsigned int numTasks;
  T* pData;
};


/**
  Creates the data stuctures that will be sent to each thread.  Divides the vertices into
  discrete chunks to be evaluated in the threads.
  @param[in] taskCount The number of individual tasks we want to divide the calculation into.
  @param[in] geomIndex The index of the input geometry we are evaluating.
*/
template <typename T>
void CreateThreadData(int taskCount, unsigned int elementCount, T* taskData, ThreadData<T>* threadData) {
  unsigned int taskLength = (elementCount + taskCount - 1) / taskCount;
  unsigned int start = 0;
  unsigned int end = taskLength;
  int lastTask = taskCount - 1;
  for(int i = 0; i < taskCount; i++) {
    if (i == lastTask) {
      end = elementCount;
    }
    threadData[i].start = start;
    threadData[i].end = end;
    threadData[i].numTasks = taskCount;
    threadData[i].pData = taskData;

    start += taskLength;
    end += taskLength;
  }
}


#endif