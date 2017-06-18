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

#ifdef __AVX__
#include <xmmintrin.h>
#include <immintrin.h>
#endif

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
  @param[in] intermediate true to get the intermediate shape.
  @return MStatus.
 */
MStatus GetShapeNode(MDagPath& path, bool intermediate=false);


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

/**
  Calculates the components necessary to create a wrap basis matrix.
  @param[in] weights The sample weights array from the wrap binding.
  @param[in] coords The barycentric coordinates of the closest point.
  @param[in] triangleVertices The vertex ids forming the triangle of the closest point.
  @param[in] points The driver point array.
  @param[in] normals The driver per-vertex normal array.
  @param[in] sampleIds The vertex ids on the driver of the current sample.
  @param[in] alignedStorage double array that is 32 byte aligned for AVX.
  @param[out] origin The origin of the coordinate system.
  @param[out] up The up vector of the coordinate system.
  @param[out] normal The normal vector of the coordinate system.
*/
void CalculateBasisComponents(const MDoubleArray& weights, const BaryCoords& coords,
                              const MIntArray& triangleVertices, const MPointArray& points,
                              const MFloatVectorArray& normals, const MIntArray& sampleIds,
                              double* alignedStorage,
                              MPoint& origin, MVector& up, MVector& normal);

/**
  Ensures that the up and normal vectors are perpendicular to each other.
  @param[in] weights The sample weights array from the wrap binding.
  @param[in] points The driver point array.
  @param[in] sampleIds The vertex ids on the driver of the current sample.
  @param[in] origin The origin of the coordinate system.
  @param[in] up The up vector of the coordinate system.
  @param[out] normal The normal vector of the coordinate system.
*/
void GetValidUp(const MDoubleArray& weights, const MPointArray& points,
                const MIntArray& sampleIds, const MPoint& origin, const MVector& normal,
                MVector& up);


template <typename T>
struct ThreadData {
  unsigned int start;
  unsigned int end;
  unsigned int numTasks;
  double* alignedStorage;
  T* pData;

#ifdef __AVX__
  ThreadData() {
    alignedStorage = (double*) _mm_malloc(4*sizeof(double), 256);
  }
  ~ThreadData() {
    _mm_free(alignedStorage);
  }
#endif
};


/**
  Creates the data stuctures that will be sent to each thread.  Divides the vertices into
  discrete chunks to be evaluated in the threads.
  @param[in] taskCount The number of individual tasks we want to divide the calculation into.
  @param[in] elementCount The number of vertices or elements to be divided up.
  @param[in] taskData The TaskData or BindData object.
  @param[out] threadData The array of ThreadData objects.  It is assumed the array is of size taskCount.
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

#ifdef __AVX__
/**
  Calculates 4 dot products at once.
  @param[in] w1 Weight vector x element.
  @param[in] w2 Weight vector y element.
  @param[in] w3 Weight vector z element.
  @param[in] w4 Weight vector w element.
  @param[in] p1 First vector.
  @param[in] p2 Second vector.
  @param[in] p3 Third vector.
  @param[in] p4 Fourth vector.
  @return A __m256d vector where each element is the corresponding p vector dot product with w.
 */
template <typename T>
__m256d Dot4(double w1, double w2, double w3, double w4,
             const T& p1, const T& p2, const T& p3, const T& p4) {
  __m256d xxx = _mm256_set_pd(p1.x, p2.x, p3.x, p4.x);
  __m256d yyy = _mm256_set_pd(p1.y, p2.y, p3.y, p4.y);
  __m256d zzz = _mm256_set_pd(p1.z, p2.z, p3.z, p4.z);
  __m256d www = _mm256_set_pd(w1, w2, w3, w4);
  __m256d xw = _mm256_mul_pd(xxx, www);
  __m256d yw = _mm256_mul_pd(yyy, www);
  __m256d zw = _mm256_mul_pd(zzz, www);
  __m256d ww = _mm256_mul_pd(www, www); // Dummy
  // low to high: xw0+xw1 yw0+yw1 xw2+xw3 yw2+yw3
  __m256d temp01 = _mm256_hadd_pd(xw, yw);   
  // low to high: zw0+zw1 ww0+ww1 zw2+zw3 ww2+ww3
  __m256d temp23 = _mm256_hadd_pd(zw, ww);
  // low to high: xw2+xw3 yw2+yw3 zw0+zw1 ww0+ww1
  __m256d swapped = _mm256_permute2f128_pd(temp01, temp23, 0x21);
  // low to high: xw0+xw1 yw0+yw1 zw2+zw3 ww2+ww3
  __m256d blended = _mm256_blend_pd(temp01, temp23, 0xC);
  // low to high: xw0+xw1+xw2+xw3 yw0+yw1+yw2+yw3 zw0+zw1+zw2+zw3 ww0+ww1+ww2+ww3
  __m256d dotproduct = _mm256_add_pd(swapped, blended);
  return dotproduct;
}
#endif

#endif
