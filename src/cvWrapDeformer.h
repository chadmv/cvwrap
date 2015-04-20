#ifndef CVWRAPDEFORMER_H
#define CVWRAPDEFORMER_H

#include <maya/MArrayDataHandle.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MDoubleArray.h>
#include <maya/MFloatArray.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MFnComponentListData.h>
#include <maya/MFnData.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MFnMesh.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MGlobal.h>
#include <maya/MIntArray.h>
#include <maya/MMatrix.h> 
#include <maya/MMatrixArray.h> 
#include <maya/MPlug.h> 
#include <maya/MPoint.h> 
#include <maya/MPointArray.h> 
#include <maya/MThreadPool.h>
#include <maya/MTransformationMatrix.h>
#include <maya/MTypeId.h> 
#include <maya/MVector.h>
#include <maya/MVectorArray.h>

#include <maya/MItGeometry.h>
#include <maya/MItMeshVertex.h>

#include <maya/MPxDeformerNode.h>


#include <map>
#include <vector>
#include "common.h"

struct TaskData {
  MMatrix driverMatrix;
  MMatrix drivenMatrix;
  MMatrix drivenInverseMatrix;
  float envelope;
  float scale;

  MIntArray membership;
  MFloatArray paintWeights;
  MPointArray points;

  MPointArray driverPoints;
  MVectorArray driverNormals;
  MFloatVectorArray driverTangents;
  MMatrixArray bindMatrices;
  std::vector<MIntArray> triangleVerts;
  std::vector<MIntArray> tangentIndices;
  std::vector<MIntArray> sampleIds;
  std::vector<BaryCoords> baryCoords;
  std::vector<MDoubleArray> sampleWeights;
};
 

struct ThreadData {
  unsigned int start;
  unsigned int end;
  unsigned int numTasks;
  TaskData* pData;
};


class CVWrap : public MPxDeformerNode {
 public:
  CVWrap();
  virtual ~CVWrap(); 
  virtual MStatus deform(MDataBlock& data, MItGeometry& iter, const MMatrix& mat,
                         unsigned int mIndex);
  virtual MStatus setDependentsDirty(const MPlug& plugBeingDirtied, MPlugArray& affectedPlugs);

  MStatus GetBindInfo(MDataBlock& data, unsigned int geomIndex);

  static void* creator();
  static MStatus initialize();

  /**
    Creates the data stuctures that will be sent to each thread.  Divides the vertices into
    discrete chunks to be evaluated in the threads.
    @param[in] taskCount The number of individual tasks we want to divide the calculation into.
    @param[in] geomIndex The index of the input geometry we are evaluating.
  */
  void CreateThreadData(int taskCount, unsigned int geomIndex);

  /**
    Distributes the ThreadData objects to the parallel threads.
    @param[in] data The user defined data.  In this case, the ThreadData array.
    @param[in] pRoot Maya's root task.
  */
  static void CreateTasks(void *data, MThreadRootTask *pRoot);
  static MThreadRetVal EvaluateWrap(void *pParam);
    
  const static char* kName;  /**< The name of the node. */
  static MObject aBindDriverGeo;
  static MObject aDriverGeo;
  static MObject aBindData;
  static MObject aSampleVerts;
  static MObject aRadius;
  static MObject aBindInfo;
  static MObject aSampleComponents;
  static MObject aSampleWeights;
  static MObject aBindMatrix;
  /** The vertex indices of the triangle containing the origin of each coordinate system. */
  static MObject aTriangleVerts;
  /** The indices of the tangents used to calculate the up vector. */
  static MObject aTangentIndices;  
  static MObject aBarycentricWeights;
  static MObject aUVSet;
  static MObject aNumTasks;
  static MObject aScale;
  static MTypeId id;

private:
  std::map<unsigned int, bool> dirty_;
  std::vector<TaskData> taskData_;  /**< Per geometry evaluation data. */
  std::vector<ThreadData*> threadData_;

};

#endif
