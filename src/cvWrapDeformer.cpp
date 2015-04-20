#include "cvWrapDeformer.h"
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnMessageAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <cassert>

MTypeId CVWrap::id(0x0011580B);

const char* CVWrap::kName = "cvWrap";
MObject CVWrap::aBindDriverGeo;
MObject CVWrap::aDriverGeo;
MObject CVWrap::aBindData;
MObject CVWrap::aSampleVerts;
MObject CVWrap::aRadius;
MObject CVWrap::aBindInfo;
MObject CVWrap::aSampleComponents;
MObject CVWrap::aSampleWeights;
MObject CVWrap::aBindMatrix;
MObject CVWrap::aTriangleVerts;
MObject CVWrap::aTangentIndices;
MObject CVWrap::aBarycentricWeights;
MObject CVWrap::aUVSet;
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

  aUVSet = tAttr.create("uvSet", "uvSet", MFnData::kString);
  addAttribute(aUVSet);

  /* Each outputGeometry needs:
  -- bindData
     | -- sampleComponents
     | -- sampleWeights
     | -- bindMatrix
     | -- triangleVerts
     | -- barycentricWeights
  */

  aSampleComponents = tAttr.create("sampleComponents", "sampleComponents", MFnData::kComponentList);
  tAttr.setArray(true);

  aSampleWeights = tAttr.create("sampleWeights", "sampleWeights", MFnData::kDoubleArray);
  tAttr.setArray(true);

  aBindMatrix = mAttr.create("bindMatrix", "bindMatrix");
  mAttr.setDefault(MMatrix::identity);
  mAttr.setArray(true);

  aTriangleVerts = nAttr.create("triangleVerts", "triangleVerts", MFnNumericData::k3Int);
  nAttr.setArray(true);

  aTangentIndices = nAttr.create("tangentIndices", "tangentIndices", MFnNumericData::k3Int);
  nAttr.setArray(true);

  aBarycentricWeights = nAttr.create("barycentricWeights", "barycentricWeights", MFnNumericData::k3Float);
  nAttr.setArray(true);

  aBindData = cAttr.create("bindData", "bindData");
  cAttr.setArray(true);
  cAttr.addChild(aSampleComponents);
  cAttr.addChild(aSampleWeights);
  cAttr.addChild(aBindMatrix);
  cAttr.addChild(aTriangleVerts);
  cAttr.addChild(aTangentIndices);
  cAttr.addChild(aBarycentricWeights);
  addAttribute(aBindData);
  attributeAffects(aSampleComponents, outputGeom);
  attributeAffects(aSampleWeights, outputGeom);
  attributeAffects(aBindMatrix, outputGeom);
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

  MGlobal::executeCommand("makePaintable -attrType multiFloat -sm deformer CVWrap weights");
    
  return MS::kSuccess;
}


CVWrap::CVWrap() {
  MThreadPool::init();
}

CVWrap::~CVWrap() {
  MThreadPool::release();
  std::vector<ThreadData*>::iterator iter;
  for (iter = threadData_.begin(); iter != threadData_.end(); ++iter) {
    delete [] *iter;
  }
  threadData_.clear();
}


void* CVWrap::creator() { return new CVWrap(); }


MStatus CVWrap::setDependentsDirty(const MPlug& plugBeingDirtied, MPlugArray& affectedPlugs) {
  // TODO: Extract the geom index from the dirty plug and set the dirty flag
  unsigned int geomIndex = 0;
  dirty_[geomIndex] = true;
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
  if (dirty_[geomIndex]) {
    dirty_[geomIndex] = false;
    status = GetBindInfo(data, geomIndex);
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
  for(unsigned int i = 0; i < numDriverPoints; i++) {
    status = fnDriver.getVertexNormal(i, false, taskData.driverNormals[i], MSpace::kWorld);
    CHECK_MSTATUS_AND_RETURN_IT(status);
  }
  // Get the driver tangents
  MString uvSet = data.inputValue(aUVSet).asString();
  status = fnDriver.getTangents(taskData.driverTangents, MSpace::kWorld, &uvSet);
  CHECK_MSTATUS_AND_RETURN_IT(status);

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
  int numTasks = data.inputValue(aNumTasks).asInt();
  if (taskData.envelope == 0.0f || numTasks <= 0) {
    return MS::kSuccess;
  }
  CreateThreadData(numTasks, geomIndex);
  MThreadPool::newParallelRegion(CreateTasks, (void *)threadData_[geomIndex]);

  status = itGeo.setAllPositions(taskData.points);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  return MS::kSuccess;
}


MStatus CVWrap::GetBindInfo(MDataBlock& data, unsigned int geomIndex) {
  MStatus status;
  MArrayDataHandle hBindDataArray = data.inputArrayValue(aBindData);
  status = hBindDataArray.jumpToElement(geomIndex);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MDataHandle hBindData = hBindDataArray.inputValue();

  MArrayDataHandle hSampleWeights = hBindData.child(aSampleWeights);
  unsigned int numVerts = hSampleWeights.elementCount();
  if (numVerts == 0) {
    // No binding information yet.
    return MS::kNotImplemented;
  }
  MArrayDataHandle hComponents = hBindData.child(aSampleComponents);
  MArrayDataHandle hBindMatrix = hBindData.child(aBindMatrix);
  MArrayDataHandle hTriangleVerts = hBindData.child(aTriangleVerts);
  MArrayDataHandle hTangentIndices = hBindData.child(aTangentIndices);
  MArrayDataHandle hBarycentricWeights = hBindData.child(aBarycentricWeights);
  hSampleWeights.jumpToArrayElement(0);
  hComponents.jumpToArrayElement(0);
  hBindMatrix.jumpToArrayElement(0);
  hTriangleVerts.jumpToArrayElement(0);
  hTangentIndices.jumpToArrayElement(0);
  hBarycentricWeights.jumpToArrayElement(0);
  MFnSingleIndexedComponent fnSingleComp;
  MFnComponentListData fnCompData;
  MFnNumericData fnNumericData;
  TaskData& taskData = taskData_[geomIndex];
  taskData.bindMatrices.setLength(numVerts);
  taskData.triangleVerts.resize(numVerts);
  taskData.tangentIndices.resize(numVerts);
  taskData.sampleIds.resize(numVerts);
  taskData.baryCoords.resize(numVerts);
  taskData.sampleWeights.resize(numVerts);
  int sampleLength = (int)taskData.bindMatrices.length();
  for (unsigned int i = 0; i < numVerts; ++i) {
    int logicalIndex = hComponents.elementIndex();
    if (logicalIndex >= sampleLength) {
      // Nurbs surfaces may be sparse so make sure we have enough space.
      taskData.bindMatrices.setLength(logicalIndex+1);
      taskData.triangleVerts.resize(logicalIndex+1);
      taskData.tangentIndices.resize(logicalIndex+1);
      taskData.sampleIds.resize(logicalIndex+1);
      taskData.baryCoords.resize(logicalIndex+1);
      taskData.sampleWeights.resize(logicalIndex+1);
    }

    // Get sample ids
    MObject oComponentData = hComponents.inputValue().data();
    status = fnCompData.setObject(oComponentData);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MIntArray& sampleIds = taskData.sampleIds[logicalIndex];
    sampleIds.clear();
    for (unsigned int j = 0; j < fnCompData.length(); j++) {
      if (fnCompData[j].hasFn(MFn::kSingleIndexedComponent)) {
        status = fnSingleComp.setObject(fnCompData[j]);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        MIntArray tempIds;
        fnSingleComp.getElements(tempIds);
        for(unsigned int k = 0; k < tempIds.length(); k++) {
          sampleIds.append(tempIds[k]);
        }
      }
    }

    // Get sample weights
    MObject oWeightData = hSampleWeights.inputValue().data();
    MFnDoubleArrayData fnDoubleData(oWeightData, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    taskData.sampleWeights[logicalIndex] = fnDoubleData.array();
    assert(taskData.sampleWeights[logicalIndex].length() == sampleIds.length());

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

    // Get tangent indices
    int3& tangentIndices = hTangentIndices.inputValue(&status).asInt3();
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MIntArray& tangentIds = taskData.tangentIndices[logicalIndex];
    tangentIds.setLength(3);
    tangentIds[0] = tangentIndices[0];
    tangentIds[1] = tangentIndices[1];
    tangentIds[2] = tangentIndices[2];

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
    hTangentIndices.next();
    hBarycentricWeights.next();
  }
  return MS::kSuccess;
}


void CVWrap::CreateThreadData(int taskCount, unsigned int geomIndex) {
  if (geomIndex >= threadData_.size()) {
    // Make sure a ThreadData objects exist for this geomIndex.
    size_t currentSize = threadData_.size();
    threadData_.resize(geomIndex+1);
    for (size_t i = currentSize; i < geomIndex+1; ++i) {
      threadData_[i] = new ThreadData[taskCount];
    }
  } else {
    // Make sure the number of ThreadData instances is correct for this geomIndex
    if (threadData_[geomIndex][0].numTasks != taskCount) {
      delete [] threadData_[geomIndex];
      threadData_[geomIndex] = new ThreadData[taskCount];
    }
  }
  TaskData& taskData = taskData_[geomIndex];
  unsigned int taskLength = (taskData.points.length() + taskCount - 1) / taskCount;
  unsigned int start = 0;
  unsigned int end = taskLength;
  int lastTask = taskCount - 1;
  for(int i = 0; i < taskCount; i++) {
    if (i == lastTask) {
      end = taskData.points.length();
    }
    threadData_[geomIndex][i].start = start;
    threadData_[geomIndex][i].end = end;
    threadData_[geomIndex][i].numTasks = taskCount;
    threadData_[geomIndex][i].pData = &taskData;

    start += taskLength;
    end += taskLength;
  }
}


void CVWrap::CreateTasks(void *data, MThreadRootTask *pRoot) {
  ThreadData* threadData = static_cast<ThreadData*>(data);

  if (threadData) {
    int numTasks = threadData[0].numTasks;
    for(int i = 0; i < numTasks; i++) {
      MThreadPool::createTask(EvaluateWrap, (void *)&threadData[i], pRoot);
    }
    MThreadPool::executeAndJoin(pRoot);
  }
}


MThreadRetVal CVWrap::EvaluateWrap(void *pParam) {
  ThreadData* pThreadData = static_cast<ThreadData*>(pParam);
  TaskData* pData = pThreadData->pData;
  // Get the data out of the struct so it is easier to work with.
  MMatrix& driverMatrix = pData->driverMatrix;
  MMatrix& drivenMatrix = pData->drivenMatrix;
  MMatrix& drivenInverseMatrix = pData->drivenInverseMatrix;
  float env = pThreadData->pData->envelope;
  float scale = pThreadData->pData->scale;
  MIntArray& membership = pData->membership;
  MFloatArray& paintWeights = pData->paintWeights;
  MPointArray& points = pData->points;
  MPointArray& driverPoints = pData->driverPoints;
  MVectorArray& driverNormals = pData->driverNormals;
  MFloatVectorArray& driverTangents = pData->driverTangents;
  MMatrixArray& bindMatrices = pData->bindMatrices;
  std::vector <MIntArray>& triangleVerts = pData->triangleVerts;
  std::vector <MIntArray>& tangentIndices = pData->tangentIndices;
  std::vector <MIntArray>& sampleIds = pData->sampleIds;
  std::vector <BaryCoords>& baryCoords = pData->baryCoords;
  std::vector <MDoubleArray>& sampleWeights = pData->sampleWeights;

  unsigned int taskStart = pThreadData->start;
  unsigned int taskEnd = pThreadData->end;

  MPoint newPt, pt;
  MMatrix scaleMatrix;
  scaleMatrix[0][0] = scale;
  scaleMatrix[1][1] = scale;
  scaleMatrix[2][2] = scale;
  for (unsigned int i = taskStart; i < taskEnd; ++i) {
    int index = membership[i];
    MIntArray& triVerts = triangleVerts[index];
    BaryCoords& coords = baryCoords[index];

    // Reconstruct origin
    MPoint origin = (MVector(driverPoints[triVerts[0]]) * coords[0]) + 
                    (MVector(driverPoints[triVerts[1]]) * coords[1]) +
                    (MVector(driverPoints[triVerts[2]]) * coords[2]);


    // Reconstruct the up vector from the tangents
    MVector up;
    up += driverTangents[tangentIndices[index][0]] * coords[0];
    up += driverTangents[tangentIndices[index][1]] * coords[1];
    up += driverTangents[tangentIndices[index][2]] * coords[2];
    up.normalize();

    // Reconstruct normal vector
    unsigned int numSamples = sampleIds[index].length();
    MVector normal;
    MDoubleArray& weights = sampleWeights[index];
    MIntArray& ids = sampleIds[index];
    unsigned int wsize = weights.length();
    unsigned int idssize = ids.length();
    for (unsigned int j = 0; j < numSamples; j++) {
      normal += driverNormals[ids[j]] * weights[j];
    }
    normal.normalize();
        
    MMatrix matrix;
    CreateMatrix(origin, normal, up, matrix);
    matrix = scaleMatrix * matrix;
    MPoint newPt = ((points[i]  * drivenMatrix) * (bindMatrices[index] * matrix)) * drivenInverseMatrix;
    points[i] = points[i] + ((newPt - points[i]) * paintWeights[i] * env);
  }
  return 0;
}

