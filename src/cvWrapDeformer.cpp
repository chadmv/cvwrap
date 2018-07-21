#include "cvWrapDeformer.h"
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMessageAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MGlobal.h>
#include <maya/MItGeometry.h>
#include <maya/MNodeMessage.h>
#include <maya/MPlugArray.h>
#include <cassert>

MTypeId CVWrap::id(0x0011580B);

const char* CVWrap::kName = "cvWrap";
MObject CVWrap::aBindDriverGeo;
MObject CVWrap::aDriverGeo;
MObject CVWrap::aBindData;
MObject CVWrap::aSampleComponents;
MObject CVWrap::aSampleWeights;
MObject CVWrap::aTriangleVerts;
MObject CVWrap::aBarycentricWeights;
MObject CVWrap::aBindMatrix;
MObject CVWrap::aNumTasks;
MObject CVWrap::aScale;

MStatus CVWrap::initialize() {
  MFnCompoundAttribute cAttr;
  MFnMatrixAttribute mAttr;
  MFnMessageAttribute meAttr;
  MFnTypedAttribute tAttr;
  MFnNumericAttribute nAttr;
  MStatus status;

  aDriverGeo = tAttr.create("driver", "driver", MFnData::kMesh);
  addAttribute(aDriverGeo);
  attributeAffects(aDriverGeo, outputGeom);

  aBindDriverGeo = meAttr.create("bindMesh", "bindMesh");
  addAttribute(aBindDriverGeo);

  /* Each outputGeometry needs:
  -- bindData
     | -- sampleComponents
     | -- sampleWeights
     | -- triangleVerts
     | -- barycentricWeights
     | -- bindMatrix
  */

  aSampleComponents = tAttr.create("sampleComponents", "sampleComponents", MFnData::kIntArray);
  tAttr.setArray(true);

  aSampleWeights = tAttr.create("sampleWeights", "sampleWeights", MFnData::kDoubleArray);
  tAttr.setArray(true);

  aTriangleVerts = nAttr.create("triangleVerts", "triangleVerts", MFnNumericData::k3Int);
  nAttr.setArray(true);

  aBarycentricWeights = nAttr.create("barycentricWeights", "barycentricWeights", MFnNumericData::k3Float);
  nAttr.setArray(true);

  aBindMatrix = mAttr.create("bindMatrix", "bindMatrix");
  mAttr.setDefault(MMatrix::identity);
  mAttr.setArray(true);

  aBindData = cAttr.create("bindData", "bindData");
  cAttr.setArray(true);
  cAttr.addChild(aSampleComponents);
  cAttr.addChild(aSampleWeights);
  cAttr.addChild(aTriangleVerts);
  cAttr.addChild(aBarycentricWeights);
  cAttr.addChild(aBindMatrix);
  addAttribute(aBindData);
  attributeAffects(aSampleComponents, outputGeom);
  attributeAffects(aSampleWeights, outputGeom);
  attributeAffects(aBindMatrix, outputGeom);
  attributeAffects(aTriangleVerts, outputGeom);
  attributeAffects(aBarycentricWeights, outputGeom);


  aScale = nAttr.create("scale", "scale", MFnNumericData::kFloat, 1.0);
  nAttr.setKeyable(true);
  addAttribute(aScale);
  attributeAffects(aScale, outputGeom);

  aNumTasks = nAttr.create("numTasks", "numTasks", MFnNumericData::kInt, 32);
  nAttr.setMin(1);
  nAttr.setMax(64);
  addAttribute(aNumTasks);
  attributeAffects(aNumTasks, outputGeom);

  MGlobal::executeCommand("makePaintable -attrType multiFloat -sm deformer cvWrap weights");
    
  return MS::kSuccess;
}


/**
  Utility method used by both the MPxDeformer and MPxGPUDeformer to pull the bind data out
  of the datablock.
  @param[in] data The node datablock.
  @param[in] geomIndex The geometry logical index.
  @param[out] taskData Bind info storage.
*/
MStatus GetBindInfo(MDataBlock& data, unsigned int geomIndex, TaskData& taskData) {
  MStatus status;
  MArrayDataHandle hBindDataArray = data.inputArrayValue(CVWrap::aBindData);
  status = hBindDataArray.jumpToElement(geomIndex);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MDataHandle hBindData = hBindDataArray.inputValue();

  MArrayDataHandle hSampleWeights = hBindData.child(CVWrap::aSampleWeights);
  unsigned int numVerts = hSampleWeights.elementCount();
  if (numVerts == 0) {
    // No binding information yet.
    return MS::kNotImplemented;
  }
  MArrayDataHandle hComponents = hBindData.child(CVWrap::aSampleComponents);
  MArrayDataHandle hBindMatrix = hBindData.child(CVWrap::aBindMatrix);
  MArrayDataHandle hTriangleVerts = hBindData.child(CVWrap::aTriangleVerts);
  MArrayDataHandle hBarycentricWeights = hBindData.child(CVWrap::aBarycentricWeights);

  hSampleWeights.jumpToArrayElement(0);
  hComponents.jumpToArrayElement(0);
  hBindMatrix.jumpToArrayElement(0);
  hTriangleVerts.jumpToArrayElement(0);
  hBarycentricWeights.jumpToArrayElement(0);

  MFnNumericData fnNumericData;
  taskData.bindMatrices.setLength(numVerts);
  taskData.sampleIds.resize(numVerts);
  taskData.sampleWeights.resize(numVerts);
  taskData.triangleVerts.resize(numVerts);
  taskData.baryCoords.resize(numVerts);

  int sampleLength = (int)taskData.bindMatrices.length();
  for (unsigned int i = 0; i < numVerts; ++i) {
    int logicalIndex = hComponents.elementIndex();
    if (logicalIndex >= sampleLength) {
      // Nurbs surfaces may be sparse so make sure we have enough space.
      taskData.bindMatrices.setLength(logicalIndex+1);
      taskData.sampleIds.resize(logicalIndex+1);
      taskData.sampleWeights.resize(logicalIndex+1);
      taskData.triangleVerts.resize(logicalIndex+1);
      taskData.baryCoords.resize(logicalIndex+1);
    }

    // Get sample ids
    MObject oIndexData = hComponents.inputValue().data();
    MFnIntArrayData fnIntData(oIndexData, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    taskData.sampleIds[logicalIndex] = fnIntData.array();

    // Get sample weights
    MObject oWeightData = hSampleWeights.inputValue().data();
    MFnDoubleArrayData fnDoubleData(oWeightData, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    taskData.sampleWeights[logicalIndex] = fnDoubleData.array();
    assert(taskData.sampleWeights[logicalIndex].length() == taskData.sampleIds[logicalIndex].length());
    
    // Get bind matrix
    taskData.bindMatrices[logicalIndex] = hBindMatrix.inputValue().asMatrix();

    // Get triangle vertex binding
    int3& verts = hTriangleVerts.inputValue(&status).asInt3();
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MIntArray& triangleVerts = taskData.triangleVerts[logicalIndex];
    triangleVerts.setLength(3);
    triangleVerts[0] = verts[0];
    triangleVerts[1] = verts[1];
    triangleVerts[2] = verts[2];

    // Get barycentric weights
    float3& baryWeights = hBarycentricWeights.inputValue(&status).asFloat3();
    CHECK_MSTATUS_AND_RETURN_IT(status);
    BaryCoords& coords = taskData.baryCoords[logicalIndex];
    coords[0] = baryWeights[0];
    coords[1] = baryWeights[1];
    coords[2] = baryWeights[2];

    hSampleWeights.next();
    hComponents.next();
    hBindMatrix.next();
    hTriangleVerts.next();
    hBarycentricWeights.next();
  }
  return MS::kSuccess;
}

MStatus GetDriverData(MDataBlock& data, TaskData& taskData) {
  MStatus status;
  // Get driver geo
  MDataHandle hDriverGeo = data.inputValue(CVWrap::aDriverGeo, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MObject oDriverGeo = hDriverGeo.asMesh();
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MFnMesh fnDriver(oDriverGeo, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  // Get the driver point positions
  status = fnDriver.getPoints(taskData.driverPoints, MSpace::kWorld);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  unsigned int numDriverPoints = taskData.driverPoints.length();
  // Get the driver normals
  taskData.driverNormals.setLength(numDriverPoints);
  fnDriver.getVertexNormals(false, taskData.driverNormals, MSpace::kWorld);
  return MS::kSuccess;
}



CVWrap::CVWrap() {
  MThreadPool::init();
  onDeleteCallbackId = 0;
}

CVWrap::~CVWrap() {
  if (onDeleteCallbackId != 0)
    MMessage::removeCallback(onDeleteCallbackId);
	
  MThreadPool::release();
  std::vector<ThreadData<TaskData>*>::iterator iter;
  for (iter = threadData_.begin(); iter != threadData_.end(); ++iter) {
    delete [] *iter;
  }
  threadData_.clear();
}


void* CVWrap::creator() { return new CVWrap(); }

void CVWrap::postConstructor()
{
  MPxDeformerNode::postConstructor();

  MStatus status = MS::kSuccess;
  MObject obj = thisMObject();
  onDeleteCallbackId = MNodeMessage::addNodeAboutToDeleteCallback(obj, aboutToDeleteCB, NULL, &status);
}

void CVWrap::aboutToDeleteCB(MObject &node, MDGModifier &modifier, void *clientData)
{
  // Find any node connected to .bindMesh and delete it with the deformer, for compatibility with wrap.
  MPlug bindPlug(node, aBindDriverGeo);
  MPlugArray bindGeometries;
  bindPlug.connectedTo(bindGeometries, true, false);
  for (unsigned int i = 0; i < bindGeometries.length(); i++) {
    MObject node = bindGeometries[i].node();
    modifier.deleteNode(node);
  }
}


MStatus CVWrap::setDependentsDirty(const MPlug& plugBeingDirtied, MPlugArray& affectedPlugs) {
  // Extract the geom index from the dirty plug and set the dirty flag so we know that we need to
  // re-read the binding data.
  if (plugBeingDirtied.isElement()) {
    MPlug parent = plugBeingDirtied.array().parent();
    if (parent == aBindData) {
      unsigned int geomIndex = parent.logicalIndex();
      dirty_[geomIndex] = true;
    }
  }
  return MS::kSuccess;
}


MStatus CVWrap::deform(MDataBlock& data, MItGeometry& itGeo, const MMatrix& localToWorldMatrix,
                       unsigned int geomIndex) {
  MStatus status;
  if (geomIndex >= taskData_.size()) {
    taskData_.resize(geomIndex+1);
  }
  TaskData& taskData = taskData_[geomIndex];
  
  // Get driver geo
  MDataHandle hDriverGeo = data.inputValue(aDriverGeo, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MObject oDriverGeo = hDriverGeo.asMesh();
  CHECK_MSTATUS_AND_RETURN_IT(status);
  if (oDriverGeo.isNull()) {
    // Without a driver mesh, we can't do anything
    return MS::kSuccess;
  }

  // Only pull bind information from the data block if it is dirty
  if (dirty_[geomIndex] || taskData.sampleIds.size() == 0) {
    dirty_[geomIndex] = false;
    status = GetBindInfo(data, geomIndex, taskData);
    if (status == MS::kNotImplemented) {
      // If no bind information is stored yet, don't do anything.
      return MS::kSuccess;
    } else if (MFAIL(status)) {
      CHECK_MSTATUS_AND_RETURN_IT(status);
    }
  }

  // Get driver geo information
  MFnMesh fnDriver(oDriverGeo, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  // Get the driver point positions
  status = fnDriver.getPoints(taskData.driverPoints, MSpace::kWorld);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  unsigned int numDriverPoints = taskData.driverPoints.length();
  // Get the driver normals
  taskData.driverNormals.setLength(numDriverPoints);
  fnDriver.getVertexNormals(false, taskData.driverNormals, MSpace::kWorld);

  // Get the deformer membership and paint weights
  unsigned int membershipCount = itGeo.count();
  taskData.membership.setLength(membershipCount);
  taskData.paintWeights.setLength(membershipCount);
  status = itGeo.allPositions(taskData.points);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  for (int i = 0; !itGeo.isDone(); itGeo.next(), i++) {
    taskData.membership[i] = itGeo.index();
    taskData.paintWeights[i] = weightValue(data, geomIndex, itGeo.index());
  }
  
  taskData.drivenMatrix = localToWorldMatrix;
  taskData.drivenInverseMatrix = localToWorldMatrix.inverse();
  
  // See if we even need to calculate anything.
  taskData.scale = data.inputValue(aScale).asFloat();
  taskData.envelope = data.inputValue(envelope).asFloat();
  int taskCount = data.inputValue(aNumTasks).asInt();
  if (taskData.envelope == 0.0f || taskCount <= 0) {
    return MS::kSuccess;
  }

  if (geomIndex >= threadData_.size()) {
    // Make sure a ThreadData objects exist for this geomIndex.
    size_t currentSize = threadData_.size();
    threadData_.resize(geomIndex+1);
    for (size_t i = currentSize; i < geomIndex+1; ++i) {
      threadData_[i] = new ThreadData<TaskData>[taskCount];
    }
  } else {
    // Make sure the number of ThreadData instances is correct for this geomIndex
    if (threadData_[geomIndex][0].numTasks != taskCount) {
      delete [] threadData_[geomIndex];
      threadData_[geomIndex] = new ThreadData<TaskData>[taskCount];
    }
  }

  CreateThreadData<TaskData>(taskCount, taskData_[geomIndex].points.length(),
                             &taskData_[geomIndex], threadData_[geomIndex]);
  MThreadPool::newParallelRegion(CreateTasks, (void *)threadData_[geomIndex]);

  status = itGeo.setAllPositions(taskData.points);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  return MS::kSuccess;
}


void CVWrap::CreateTasks(void *data, MThreadRootTask *pRoot) {
  ThreadData<TaskData>* threadData = static_cast<ThreadData<TaskData>*>(data);

  if (threadData) {
    int numTasks = threadData[0].numTasks;
    for(int i = 0; i < numTasks; i++) {
      MThreadPool::createTask(EvaluateWrap, (void *)&threadData[i], pRoot);
    }
    MThreadPool::executeAndJoin(pRoot);
  }
}


MThreadRetVal CVWrap::EvaluateWrap(void *pParam) {
  ThreadData<TaskData>* pThreadData = static_cast<ThreadData<TaskData>*>(pParam);
  double*& alignedStorage = pThreadData->alignedStorage;
  TaskData* pData = pThreadData->pData;
  // Get the data out of the struct so it is easier to work with.
  MMatrix& drivenMatrix = pData->drivenMatrix;
  MMatrix& drivenInverseMatrix = pData->drivenInverseMatrix;
  float env = pThreadData->pData->envelope;
  float scale = pThreadData->pData->scale;
  MIntArray& membership = pData->membership;
  MFloatArray& paintWeights = pData->paintWeights;
  MPointArray& points = pData->points;
  MPointArray& driverPoints = pData->driverPoints;
  MFloatVectorArray& driverNormals = pData->driverNormals;
  MMatrixArray& bindMatrices = pData->bindMatrices;
  std::vector <MIntArray>& sampleIds = pData->sampleIds;
  std::vector <MDoubleArray>& sampleWeights = pData->sampleWeights;
  std::vector <MIntArray>& triangleVerts = pData->triangleVerts;
  std::vector <BaryCoords>& baryCoords = pData->baryCoords;

  unsigned int taskStart = pThreadData->start;
  unsigned int taskEnd = pThreadData->end;

  MPoint newPt;
  MMatrix scaleMatrix, matrix;
  scaleMatrix[0][0] = scale;
  scaleMatrix[1][1] = scale;
  scaleMatrix[2][2] = scale;
  for (unsigned int i = taskStart; i < taskEnd; ++i) {
    if (i >= points.length()) {
      break;
    }
    int index = membership[i];

    MPoint origin;
    MVector normal, up;
    CalculateBasisComponents(sampleWeights[index], baryCoords[index], triangleVerts[index],
                             driverPoints, driverNormals, sampleIds[index], alignedStorage,
                             origin, up, normal);

    CreateMatrix(origin, normal, up, matrix);
    matrix = scaleMatrix * matrix;
    MPoint newPt = ((points[i]  * drivenMatrix) * (bindMatrices[index] * matrix)) * drivenInverseMatrix;
    points[i] = points[i] + ((newPt - points[i]) * paintWeights[i] * env);
  }
  return 0;
}


#if MAYA_API_VERSION >= 201600
MString CVWrapGPU::pluginLoadPath;

#if MAYA_API_VERSION >= 201650
cl_command_queue (*getMayaDefaultOpenCLCommandQueue)() = MOpenCLInfo::getMayaDefaultOpenCLCommandQueue;
#else
cl_command_queue (*getMayaDefaultOpenCLCommandQueue)() = MOpenCLInfo::getOpenCLCommandQueue;
#endif
/**
  Convenience function to copy array data to the gpu.
*/
cl_int EnqueueBuffer(MAutoCLMem& mclMem, size_t bufferSize, void* data) {
  cl_int err = CL_SUCCESS;
  if (!mclMem.get())	{
    // The buffer doesn't exist yet so create it and copy the data over.
		mclMem.attach(clCreateBuffer(MOpenCLInfo::getOpenCLContext(),
                                        CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY,
                                        bufferSize, data, &err));
	}	else {
		// The buffer already exists so just copy the data over.
		err = clEnqueueWriteBuffer(getMayaDefaultOpenCLCommandQueue(),
                               mclMem.get(), CL_TRUE, 0, bufferSize,
                               data, 0, NULL, NULL);
	}
  return err;
}

MGPUDeformerRegistrationInfo* CVWrapGPU::GetGPUDeformerInfo() {
  static CVWrapGPUDeformerInfo wrapInfo;
  return &wrapInfo;
}

CVWrapGPU::CVWrapGPU() {
  // Remember the ctor must be fast.  No heavy work should be done here.
  // Maya may allocate one of these and then never use it.
}

CVWrapGPU::~CVWrapGPU() {
  terminate();
}

#if MAYA_API_VERSION <= 201700
MPxGPUDeformer::DeformerStatus CVWrapGPU::evaluate(MDataBlock& block,
                                                   const MEvaluationNode& evaluationNode,
                                                   const MPlug& plug,
                                                   unsigned int numElements,
                                                   const MAutoCLMem inputBuffer,
                                                   const MAutoCLEvent inputEvent,
                                                   MAutoCLMem outputBuffer,
                                                   MAutoCLEvent& outputEvent) {
#else
MPxGPUDeformer::DeformerStatus  CVWrapGPU::evaluate(MDataBlock& block,
													const MEvaluationNode& evaluationNode,
													const MPlug& plug,
													const MGPUDeformerData& inputData,
													MGPUDeformerData& outputData) {
	// get the input GPU data and event
	MGPUDeformerBuffer inputDeformerBuffer = inputData.getBuffer(sPositionsName());
	const MAutoCLMem inputBuffer = inputDeformerBuffer.buffer();
	unsigned int numElements = inputDeformerBuffer.elementCount();
	const MAutoCLEvent inputEvent = inputDeformerBuffer.bufferReadyEvent();

	// create the output buffer
	MGPUDeformerBuffer outputDeformerBuffer = createOutputBuffer(inputDeformerBuffer);
	MAutoCLEvent outputEvent;
	MAutoCLMem outputBuffer = outputDeformerBuffer.buffer();
#endif

  MStatus status;
  numElements_ = numElements;
  // Copy all necessary data to the gpu.
  status = EnqueueBindData(block, evaluationNode, plug);
  CHECK_MSTATUS(status);
  status = EnqueueDriverData(block, evaluationNode, plug);
  CHECK_MSTATUS(status);
  status = EnqueuePaintMapData(block, evaluationNode, numElements, plug);
  CHECK_MSTATUS(status);

  if (!kernel_.get())  {
    // Load the OpenCL kernel if we haven't yet.
    MString openCLKernelFile(pluginLoadPath);
#if MAYA_API_VERSION > 201700
	openCLKernelFile += "/cvwrap.cl";
#else
    openCLKernelFile += "/cvwrap_pre2018.cl";
#endif
    kernel_ = MOpenCLInfo::getOpenCLKernel(openCLKernelFile, "cvwrap");
    if (kernel_.isNull())  {
      std::cerr << "Could not compile kernel " << openCLKernelFile.asChar() << "\n";
      return MPxGPUDeformer::kDeformerFailure;
    }
  }
  float envelope = block.inputValue(MPxDeformerNode::envelope, &status).asFloat();
  CHECK_MSTATUS(status);
  cl_int err = CL_SUCCESS;
  
  // Set all of our kernel parameters.  Input buffer and output buffer may be changing every frame
  // so always set them.
  unsigned int parameterId = 0;
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_mem), (void*)outputBuffer.getReadOnlyRef());
  MOpenCLInfo::checkCLErrorStatus(err);
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_mem), (void*)inputBuffer.getReadOnlyRef());
  MOpenCLInfo::checkCLErrorStatus(err);
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_mem), (void*)driverPoints_.getReadOnlyRef());
  MOpenCLInfo::checkCLErrorStatus(err);
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_mem), (void*)driverNormals_.getReadOnlyRef());
  MOpenCLInfo::checkCLErrorStatus(err);
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_mem), (void*)paintWeights_.getReadOnlyRef());
  MOpenCLInfo::checkCLErrorStatus(err);
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_mem), (void*)sampleCounts_.getReadOnlyRef());
  MOpenCLInfo::checkCLErrorStatus(err);
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_mem), (void*)sampleOffsets_.getReadOnlyRef());
  MOpenCLInfo::checkCLErrorStatus(err);
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_mem), (void*)sampleIds_.getReadOnlyRef());
  MOpenCLInfo::checkCLErrorStatus(err);
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_mem), (void*)sampleWeights_.getReadOnlyRef());
  MOpenCLInfo::checkCLErrorStatus(err);
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_mem), (void*)triangleVerts_.getReadOnlyRef());
  MOpenCLInfo::checkCLErrorStatus(err);
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_mem), (void*)baryCoords_.getReadOnlyRef());
  MOpenCLInfo::checkCLErrorStatus(err);
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_mem), (void*)bindMatrices_.getReadOnlyRef());
  MOpenCLInfo::checkCLErrorStatus(err);
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_mem), (void*)drivenMatrices_.getReadOnlyRef());
  MOpenCLInfo::checkCLErrorStatus(err);
#if MAYA_API_VERSION > 201700
  // get the world space and inverse world space matrix mem handles
  MGPUDeformerBuffer inputWorldSpaceMatrixDeformerBuffer = inputData.getBuffer(sGeometryMatrixName());
  const MAutoCLMem deformerWorldSpaceMatrix = inputWorldSpaceMatrixDeformerBuffer.buffer();
  MGPUDeformerBuffer inputInvWorldSpaceMatrixDeformerBuffer = inputData.getBuffer(sInverseGeometryMatrixName());
  const MAutoCLMem deformerInvWorldSpaceMatrix = inputInvWorldSpaceMatrixDeformerBuffer.buffer();
  // Note: these matrices are in row major order
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_mem), (void*)deformerWorldSpaceMatrix.getReadOnlyRef());
  MOpenCLInfo::checkCLErrorStatus(err);
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_mem), (void*)deformerInvWorldSpaceMatrix.getReadOnlyRef());
  MOpenCLInfo::checkCLErrorStatus(err);
#endif
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_float), (void*)&envelope);
  MOpenCLInfo::checkCLErrorStatus(err);
  err = clSetKernelArg(kernel_.get(), parameterId++, sizeof(cl_uint), (void*)&numElements_);
  MOpenCLInfo::checkCLErrorStatus(err);

  // Figure out a good work group size for our kernel.
  size_t workGroupSize;
  size_t retSize;
  err = clGetKernelWorkGroupInfo(
    kernel_.get(),
    MOpenCLInfo::getOpenCLDeviceId(),
    CL_KERNEL_WORK_GROUP_SIZE,
    sizeof(size_t),
    &workGroupSize,
    &retSize);
  MOpenCLInfo::checkCLErrorStatus(err);

  size_t localWorkSize = 256;
  if (retSize > 0) {
    localWorkSize = workGroupSize;
  }
  // global work size must be a multiple of localWorkSize
  size_t globalWorkSize = (localWorkSize - numElements_ % localWorkSize) + numElements_;

  // set up our input events.  The input event could be NULL, in that case we need to pass
  // slightly different parameters into clEnqueueNDRangeKernel
  unsigned int numInputEvents = 0;
  if (inputEvent.get()) {
    numInputEvents = 1;
  }

  // run the kernel
  err = clEnqueueNDRangeKernel(
    getMayaDefaultOpenCLCommandQueue(),
    kernel_.get(),
    1,
    NULL,
    &globalWorkSize,
    &localWorkSize,
    numInputEvents,
    numInputEvents ? inputEvent.getReadOnlyRef() : 0,
    outputEvent.getReferenceForAssignment() );
  MOpenCLInfo::checkCLErrorStatus(err);

#if MAYA_API_VERSION > 201700
  // set the buffer into the output data
  outputDeformerBuffer.setBufferReadyEvent(outputEvent);
  outputData.setBuffer(outputDeformerBuffer);
#endif

  return MPxGPUDeformer::kDeformerSuccess;
}

MStatus CVWrapGPU::EnqueueBindData(MDataBlock& data, const MEvaluationNode& evaluationNode,
                                   const MPlug& plug) {
  MStatus status;
	if ((bindMatrices_.get() && (
        !evaluationNode.dirtyPlugExists(CVWrap::aBindData, &status) &&
        !evaluationNode.dirtyPlugExists(CVWrap::aSampleComponents, &status) &&
        !evaluationNode.dirtyPlugExists(CVWrap::aSampleWeights, &status) &&
        !evaluationNode.dirtyPlugExists(CVWrap::aTriangleVerts, &status) &&
        !evaluationNode.dirtyPlugExists(CVWrap::aBarycentricWeights, &status) &&
        !evaluationNode.dirtyPlugExists(CVWrap::aBindMatrix, &status)
      )) || !status) {
    // No bind data has changed, nothing to do.
    return MS::kSuccess;
  }

  TaskData taskData;
  unsigned int geomIndex = plug.logicalIndex();
  status = GetBindInfo(data, geomIndex, taskData);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  // Flatten out bind matrices to float array
  size_t arraySize = taskData.bindMatrices.length() * 16;
	float* bindMatrices = new float[arraySize];
  for(unsigned int i = 0, idx = 0; i < taskData.bindMatrices.length(); ++i) {
    for(unsigned int row = 0; row < 4; row++) {
		  for(unsigned int column = 0; column < 4; column++) {
			  bindMatrices[idx++] = (float)taskData.bindMatrices[i](row, column);
      }
		}
	}
  cl_int err = EnqueueBuffer(bindMatrices_, arraySize * sizeof(float), (void*)bindMatrices);
  delete [] bindMatrices;

  // Store samples per vertex
  arraySize = taskData.sampleIds.size();
  int* samplesPerVertex = new int[arraySize];
  int* sampleOffsets = new int[arraySize];
  int totalSamples = 0;
  for(size_t i = 0; i < taskData.sampleIds.size(); ++i) {
    samplesPerVertex[i] = (int)taskData.sampleIds[i].length();
    sampleOffsets[i] = totalSamples;
    totalSamples += samplesPerVertex[i];
  }
  err = EnqueueBuffer(sampleCounts_, arraySize * sizeof(int), (void*)samplesPerVertex);
  err = EnqueueBuffer(sampleOffsets_, arraySize * sizeof(int), (void*)sampleOffsets);
  delete [] samplesPerVertex;
  delete [] sampleOffsets;

  // Store sampleIds and sampleWeights
  int* sampleIds = new int[totalSamples];
  float* sampleWeights = new float[totalSamples];
  int iter = 0;
  for(size_t i = 0; i < taskData.sampleIds.size(); ++i) {
    for(unsigned int j = 0; j < taskData.sampleIds[i].length(); ++j) {
      sampleIds[iter] = taskData.sampleIds[i][j];
      sampleWeights[iter] = (float)taskData.sampleWeights[i][j];
      iter++;
    }
  }
  err = EnqueueBuffer(sampleIds_, totalSamples * sizeof(int), (void*)sampleIds);
  err = EnqueueBuffer(sampleWeights_, totalSamples * sizeof(float), (void*)sampleWeights);
  delete [] sampleIds;
  delete [] sampleWeights;

  // Store triangle verts and bary coords
  arraySize = taskData.triangleVerts.size() * 3;
  int* triangleVerts = new int[arraySize];
  float* baryCoords = new float[arraySize];
  iter = 0;
  for(size_t i = 0; i < taskData.triangleVerts.size(); ++i) {
    for(unsigned int j = 0; j < 3; ++j) {
      triangleVerts[iter] = taskData.triangleVerts[i][j];
      baryCoords[iter] = (float)taskData.baryCoords[i][j];
      iter++;
    }
  }
  err = EnqueueBuffer(triangleVerts_, arraySize * sizeof(int), (void*)triangleVerts);
  err = EnqueueBuffer(baryCoords_, arraySize * sizeof(float), (void*)baryCoords);
  delete [] triangleVerts;
  delete [] baryCoords;
  return MS::kSuccess;
}


MStatus CVWrapGPU::EnqueueDriverData(MDataBlock& data, const MEvaluationNode& evaluationNode, const MPlug& plug) {
  MStatus status;
  TaskData taskData;
  status = GetDriverData(data, taskData);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  cl_int err = CL_SUCCESS;
  // Store world space driver points and normals into float arrays.
  // Reuse the same array for points and normals so we're not dynamically allocating double
  // the memory.
  unsigned int pointCount = taskData.driverPoints.length();
  float* driverData = new float[pointCount * 3];

  // Store the driver points on the gpu.
  for (unsigned int i = 0, iter = 0; i < pointCount; ++i) {
    driverData[iter++] = (float)taskData.driverPoints[i].x;
    driverData[iter++] = (float)taskData.driverPoints[i].y;
    driverData[iter++] = (float)taskData.driverPoints[i].z;
  }
  err = EnqueueBuffer(driverPoints_, pointCount * 3 * sizeof(float), (void*)driverData);

  // Store the driver normals on the gpu.
  for (unsigned int i = 0, iter = 0; i < pointCount; ++i) {
    driverData[iter++] = taskData.driverNormals[i].x;
    driverData[iter++] = taskData.driverNormals[i].y;
    driverData[iter++] = taskData.driverNormals[i].z;
  }
  err = EnqueueBuffer(driverNormals_, pointCount * 3 * sizeof(float), (void*)driverData);
	delete [] driverData;

  int idx = 0;
#if MAYA_API_VERSION <= 201700
  // Store the driven matrices on the gpu.
  MArrayDataHandle hInputs = data.inputValue(CVWrap::input, &status);
  unsigned int geomIndex = plug.logicalIndex();
  status = hInputs.jumpToElement(geomIndex);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MDataHandle hInput = hInputs.inputValue(&status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MDataHandle hGeom = hInput.child(CVWrap::inputGeom);
  MMatrix localToWorldMatrix = hGeom.geometryTransformMatrix();
  MMatrix worldToLocalMatrix = localToWorldMatrix.inverse();
  float drivenMatrices[48]; // 0-15: localToWorld, 16-31: worldToLocal, 32-47: scale

  // Store in column order so we can dot in the cl kernel.
  for(unsigned int column = 0; column < 4; column++) {
    for(unsigned int row = 0; row < 4; row++) {
			drivenMatrices[idx++] = (float)localToWorldMatrix(row, column);
    }
	}
  for(unsigned int column = 0; column < 4; column++) {
    for(unsigned int row = 0; row < 4; row++) {
			drivenMatrices[idx++] = (float)worldToLocalMatrix(row, column);
    }
	}
#else
	float drivenMatrices[16]; // 0-15: scale
#endif
  // Scale matrix is stored row major
  float scale = data.inputValue(CVWrap::aScale, &status).asFloat();
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MMatrix scaleMatrix;
  scaleMatrix[0][0] = scale;
  scaleMatrix[1][1] = scale;
  scaleMatrix[2][2] = scale;
  for(unsigned int row = 0; row < 4; row++) {
    for(unsigned int column = 0; column < 4; column++) {
			drivenMatrices[idx++] = (float)scaleMatrix(row, column);
    }
	}
#if MAYA_API_VERSION <= 201700
  err = EnqueueBuffer(drivenMatrices_, 48 * sizeof(float), (void*)drivenMatrices);
#else
  err = EnqueueBuffer(drivenMatrices_, 16 * sizeof(float), (void*)drivenMatrices);
#endif
  return MS::kSuccess;
}


MStatus CVWrapGPU::EnqueuePaintMapData(MDataBlock& data,
                                       const MEvaluationNode& evaluationNode,
                                       unsigned int numElements,
                                       const MPlug& plug) {
  MStatus status;
  if ((paintWeights_.get() &&
       !evaluationNode.dirtyPlugExists(MPxDeformerNode::weightList, &status)) || !status) {
    // The paint weights are not dirty so no need to get them.
		return MS::kSuccess;
	}

  cl_int err = CL_SUCCESS;

  // Store the paint weights on the gpu.
  // Since we can't call MPxDeformerNode::weightValue, get the paint weights from the data block.
  float* paintWeights = new float[numElements];
  MArrayDataHandle weightList = data.outputArrayValue(MPxDeformerNode::weightList, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  unsigned int geomIndex = plug.logicalIndex();
  status = weightList.jumpToElement(geomIndex);
  // it is possible that the jumpToElement fails.  In that case all weights are 1.
  if (!status) {  
    for(unsigned int i = 0; i < numElements; i++) {
      paintWeights[i] = 1.0f;
    }
  } else {
    // Initialize all weights to 1.0f
    for(unsigned int i = 0; i < numElements; i++) {
      paintWeights[i] = 1.0f;
    }
    MDataHandle weightsStructure = weightList.inputValue(&status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MArrayDataHandle weights = weightsStructure.child(MPxDeformerNode::weights);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    // Gather all the non-zero weights
    unsigned int numWeights = weights.elementCount(&status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    for (unsigned int i = 0; i < numWeights; i++, weights.next()) {
      unsigned int weightsElementIndex = weights.elementIndex(&status);
      MDataHandle value = weights.inputValue(&status);
      // BUG: The weightsElementIndex may be sparse for nurbs surfaces so this would be incorrect
      paintWeights[weightsElementIndex] = value.asFloat();
    }
  }
  err = EnqueueBuffer(paintWeights_, numElements * sizeof(float), (void*)paintWeights);
  delete [] paintWeights;
  return MS::kSuccess;
}


void CVWrapGPU::terminate() {
  driverPoints_.reset();
  driverNormals_.reset();
  paintWeights_.reset();
  bindMatrices_.reset();
  sampleCounts_.reset();
  sampleIds_.reset();
  sampleWeights_.reset();
  triangleVerts_.reset();
  baryCoords_.reset();
  drivenMatrices_.reset();
	MOpenCLInfo::releaseOpenCLKernel(kernel_);
	kernel_.reset();
}

#endif

