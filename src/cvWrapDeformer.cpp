#include "cvWrapDeformer.h"
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnMessageAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MThreadPool.h>

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

  aBarycentricWeights = nAttr.create("barycentricWeights", "barycentricWeights", MFnNumericData::k3Float);
  nAttr.setArray(true);

  aBindData = cAttr.create("bindData", "bindData");
  cAttr.setArray(true);
  cAttr.addChild(aSampleComponents);
  cAttr.addChild(aSampleWeights);
  cAttr.addChild(aBindMatrix);
  cAttr.addChild(aTriangleVerts);
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

  m_triangleVerts.clear();
  m_barycentricWeights.clear();
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
  MArrayDataHandle hBarycentricWeights = hBindData.child(aBarycentricWeights);
  hSampleWeights.jumpToElement(0);
  hComponents.jumpToElement(0);
  hBindMatrix.jumpToElement(0);
  hTriangleVerts.jumpToElement(0);
  hBarycentricWeights.jumpToElement(0);
  MFnSingleIndexedComponent fnSingleComp;
  MFnComponentListData fnCompData;
  MFnNumericData fnNumericData;
  TaskData& taskData = taskData_[geomIndex];
  taskData.bindMatrices.setLength(numVerts);
  taskData.triangleVerts.resize(numVerts);
  taskData.sampleIds.resize(numVerts);
  taskData.baryCoords.resize(numVerts);
  taskData.sampleWeights.resize(numVerts);
  for (unsigned int i = 0; i < numVerts; ++i) {
    int logicalIndex = hComponents.elementIndex();
    if (logicalIndex >= taskData.bindMatrices.length()) {
      // Nurbs surfaces may be sparse so make sure we have enough space.
      taskData.bindMatrices.setLength(logicalIndex+1);
      taskData.triangleVerts.resize(logicalIndex+1);
      taskData.sampleIds.resize(logicalIndex+1);
      taskData.baryCoords.resize(logicalIndex+1);
      taskData.sampleWeights.resize(logicalIndex+1);
    }

    // Get sample ids
    MObject oComponentData = hComponents.inputValue().data();
    status = fnCompData.setObject(oComponentData);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MIntArray& sampleIds = taskData.sampleIds[logicalIndex];
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


void CVWrap::CreateThreadData(int numTasks, unsigned int geomIndex) {
  if (geomIndex >= threadData_.size()) {
    // Make sure a ThreadData object exists for this geomIndex.
    size_t currentSize = threadData_.size();
    threadData_.resize(geomIndex+1);
    for (size_t i = currentSize; i < geomIndex+1; ++i) {
      threadData_[i] = new ThreadData[numTasks];
    }
  } else {
    // Make sure the number of ThreadData instances is correct for this geomIndex
    if (threadData_[geomIndex][0].numTasks != numTasks) {
      delete [] threadData_[geomIndex];
      threadData_[geomIndex] = new ThreadData[numTasks];
    }
  }
  TaskData& taskData = taskData_[geomIndex];
  unsigned int taskLength = (taskData.points.length() + numTasks - 1) / numTasks;
  unsigned int start = 0;
  unsigned int end = taskLength;
  int lastTask = numTasks - 1;
  for(int i = 0; i < numTasks; i++) {
    if (i == lastTask) {
      end = taskData.points.length();
    }
    threadData_[geomIndex][i].start = start;
    threadData_[geomIndex][i].end = end;
    threadData_[geomIndex][i].numTasks = numTasks;
    threadData_[geomIndex][i].pData = &taskData;

    start += taskLength;
    end += taskLength;
  }
}


void CVWrap::CreateTasks(void *pData, MThreadRootTask *pRoot) {
  ThreadData* threadData = static_cast<ThreadData*>(pData);

  if (threadData) {
    int numTasks = threadData[0].numTasks;
    for(int i = 0; i < numTasks; i++) {
      MThreadPool::createTask(threadEvaluate, (void *)&threadData[i], pRoot);
    }
    MThreadPool::executeAndJoin(pRoot);
  }
}


MThreadRetVal CVWrap::threadEvaluate(void *pParam)
{
    PTHREADDATA pThreadData = (PTHREADDATA)(pParam);
    PTASKDATA pData = pThreadData->pData;
    float env = pThreadData->pData->envelope;
    MIntArray& membership = pData->membership;
    MFloatArray& weights = pData->weights;
    MPointArray& points = pData->points;
    unsigned int numDeformVerts = points.length();

    MPointArray& driverPoints = pData->driverPoints;
    MVectorArray& driverNormals = pData->driverNormals;
    std::vector < std::vector <Sortable> >& sortedWeights = pData->sortedWeights;
    MMatrixArray& bindMatrices = pData->bindMatrices;
    std::vector <MIntArray>& triangleVerts = pData->triangleVerts;
    std::vector <MFloatArray>& barycentricWeights = pData->barycentricWeights;
    MMatrix& driverMatrix = pData->driverMatrix;
    MMatrix& drivenMatrix = pData->drivenMatrix;
    MMatrix& drivenInverseMatrix = pData->drivenInverseMatrix;
    float scale = pThreadData->pData->scale;

    unsigned int taskStart = pThreadData->start;
    unsigned int taskEnd = pThreadData->end;
    int index;

    unsigned int numSamples;
    MPoint origin, newPt, pt;
    MVector normal, up;
    MMatrix matrix;
  for (unsigned int i = taskStart; i < taskEnd; i++)
    {
        if (i >= numDeformVerts)
        {
            continue;
        }
        index = membership[i];
        numSamples = sortedWeights[index].size();

        // Reconstruct origin
        origin = (MVector(driverPoints[triangleVerts[index][0]]) * barycentricWeights[index][0]) + 
            (MVector(driverPoints[triangleVerts[index][1]]) * barycentricWeights[index][1]) +
            (MVector(driverPoints[triangleVerts[index][2]]) * barycentricWeights[index][2]);


        // Reconstruct normal and up vector
        normal = MVector::zero;
        up = MVector::zero;
        for (unsigned int j = 0; j < numSamples; j++)
        {
            normal += driverNormals[sortedWeights[index][j].index] * sortedWeights[index][j].normalizedWeight;
            up += (driverPoints[sortedWeights[index][j].index] - origin) * sortedWeights[index][j].normalizedWeight;
        }
        normal.normalize();
        normal *= scale;

        // If the up and normal or parallel or the up has 0 length, take out some influence.
        if (up * normal == 1.0 || up.length() < 0.0001)
        {
            for (unsigned int j = 0; j < numSamples; j++)
            {
                if (up * normal != 1.0 && up.length() > 0.0001)
                {
                    break;
                }
                up -= (driverPoints[sortedWeights[index][j].index] - origin) * sortedWeights[index][j].normalizedWeight;
            }
        }   

        up.normalize();
        
        createMatrix(matrix, origin, normal, up);
        
        pt = points[i];
        newPt = ((pt  * drivenMatrix) * (bindMatrices[index] * matrix)) * drivenInverseMatrix;
        points[i] = pt + ((newPt - pt) * weights[i] * env);
    }
    return 0;
}


void CVWrap::createMatrix(MMatrix& matrix, MPoint& origin, MVector& normal, MVector& up)
{
    MVector x = normal ^ up;
    up = normal ^ x;
    matrix[0][0] = x.x;         matrix[0][1] = x.y;         matrix[0][2] = x.z;      matrix[0][3] = 0.0;
    matrix[1][0] = normal.x;    matrix[1][1] = normal.y;    matrix[1][2] = normal.z; matrix[1][3] = 0.0;
    matrix[2][0] = up.x;        matrix[2][1] = up.y;        matrix[2][2] = up.z;     matrix[2][3] = 0.0;
    matrix[3][0] = origin.x;    matrix[3][1] = origin.y;    matrix[3][2] = origin.z; matrix[3][3] = 1.0;
}





