#ifndef CVWRAPDEFORMER_H
#define CVWRAPDEFORMER_H

#include <maya/MDGModifier.h>
#include <maya/MFloatArray.h>
#include <maya/MIntArray.h>
#include <maya/MMatrix.h> 
#include <maya/MMatrixArray.h> 
#include <maya/MMessage.h>
#include <maya/MPoint.h> 
#include <maya/MThreadPool.h>
#include <maya/MPxDeformerNode.h>

#if MAYA_API_VERSION >= 201600
#include <maya/MPxGPUDeformer.h>
#include <maya/MGPUDeformerRegistry.h>
#include <maya/MOpenCLInfo.h>
#include <clew/clew_cl.h>
#endif

#include <map>
#include <vector>
#include "common.h"

struct TaskData {
  MMatrix drivenMatrix;
  MMatrix drivenInverseMatrix;
  float envelope;
  float scale;

  MIntArray membership;
  MFloatArray paintWeights;
  MPointArray points;

  MPointArray driverPoints;
  MFloatVectorArray driverNormals;
  MMatrixArray bindMatrices;
  std::vector<MIntArray> sampleIds;
  std::vector<MDoubleArray> sampleWeights;
  std::vector<MIntArray> triangleVerts;
  std::vector<BaryCoords> baryCoords;
};
 

class CVWrap : public MPxDeformerNode {
 public:
  CVWrap();
  virtual ~CVWrap(); 
  virtual void postConstructor();
  virtual MStatus deform(MDataBlock& data, MItGeometry& iter, const MMatrix& mat,
                         unsigned int mIndex);
  virtual MStatus setDependentsDirty(const MPlug& plugBeingDirtied, MPlugArray& affectedPlugs);

  static void* creator();
  static MStatus initialize();


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
  static MObject aSampleComponents;
  static MObject aSampleWeights;
    /** The vertex indices of the triangle containing the origin of each coordinate system. */
  static MObject aTriangleVerts;
  /** The indices of the tangents used to calculate the up vector. */
  static MObject aBarycentricWeights;

  static MObject aBindMatrix;
  static MObject aNumTasks;
  static MObject aScale;
  static MTypeId id;

private:
  static void aboutToDeleteCB(MObject &node, MDGModifier &modifier, void *clientData);

  std::map<unsigned int, bool> dirty_;
  std::vector<TaskData> taskData_;  /**< Per geometry evaluation data. */
  std::vector<ThreadData<TaskData>*> threadData_;
  MCallbackId onDeleteCallbackId;
};



#if MAYA_API_VERSION >= 201600
// the GPU override implementation of the offsetNode
// 

class CVWrapGPU : public MPxGPUDeformer {
 public:
	// Virtual methods from MPxGPUDeformer
	CVWrapGPU();
	virtual ~CVWrapGPU();

#if MAYA_API_VERSION <= 201700
	virtual MPxGPUDeformer::DeformerStatus evaluate(MDataBlock& block, const MEvaluationNode&,
                                                  const MPlug& plug, unsigned int numElements,
                                                  const MAutoCLMem, const MAutoCLEvent,
                                                  MAutoCLMem, MAutoCLEvent&);
#else
	virtual MPxGPUDeformer::DeformerStatus evaluate(MDataBlock& block, const MEvaluationNode& evaluationNode,
													const MPlug& plug, const MGPUDeformerData& inputData,
													MGPUDeformerData& outputData);
#endif
	virtual void terminate();

	static MGPUDeformerRegistrationInfo* GetGPUDeformerInfo();
	static bool ValidateNode(MDataBlock& block, const MEvaluationNode&, const MPlug& plug, MStringArray* messages);
  /**< The path of where the plug-in is loaded from.  Used to find the cl kernel. */
  static MString pluginLoadPath;

private:
	// helper methods
	MStatus EnqueueBindData(MDataBlock& data, const MEvaluationNode& evaluationNode, const MPlug& plug);
	MStatus EnqueueDriverData(MDataBlock& data, const MEvaluationNode& evaluationNode, const MPlug& plug);
	MStatus EnqueuePaintMapData(MDataBlock& data, const MEvaluationNode& evaluationNode, unsigned int numElements, const MPlug& plug);

	// Storage for data on the GPU
	MAutoCLMem driverPoints_;
	MAutoCLMem driverNormals_;
	MAutoCLMem paintWeights_;
	MAutoCLMem bindMatrices_;
	MAutoCLMem sampleCounts_;
	MAutoCLMem sampleOffsets_;
	MAutoCLMem sampleIds_;
	MAutoCLMem sampleWeights_;
	MAutoCLMem triangleVerts_;
	MAutoCLMem baryCoords_;
	MAutoCLMem drivenMatrices_;

	unsigned int numElements_;

	// Kernel
	MAutoCLKernel kernel_;
};


/**
  The 
*/
class CVWrapGPUDeformerInfo : public MGPUDeformerRegistrationInfo {
 public:
	CVWrapGPUDeformerInfo(){}
	virtual ~CVWrapGPUDeformerInfo(){}

	virtual MPxGPUDeformer* createGPUDeformer()	{
		return new CVWrapGPU();
	}
	


#if MAYA_API_VERSION >= 201650
	virtual bool validateNodeInGraph(MDataBlock& block, const MEvaluationNode& evaluationNode,
                                   const MPlug& plug, MStringArray* messages)	{
		return true;
	}

	virtual bool validateNodeValues(MDataBlock& block, const MEvaluationNode& evaluationNode,
                                  const MPlug& plug, MStringArray* messages) {
		return true;
	}
#else
  virtual bool validateNode(MDataBlock& block, const MEvaluationNode& evaluationNode,
                            const MPlug& plug, MStringArray* messages) {
		return true;
	}
#endif
};

#endif // End Maya 2016

#endif
