#include "cvWrapCmd.h"
#include "cvWrapDeformer.h"
#include "bindingio.h"

#include <maya/MArgDatabase.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MFnMatrixData.h>
#include <maya/MFnMesh.h>
#include <maya/MGlobal.h>
#include <maya/MItDependencyGraph.h>
#include <maya/MItGeometry.h>
#include <maya/MItSelectionList.h>
#include <maya/MMeshIntersector.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MFnWeightGeometryFilter.h>
#include <maya/MSyntax.h>
#include <algorithm>
#include <cassert>
#include <utility>

#define PROGRESS_STEP 100
#define TASK_COUNT 32

/**
  A version number used to support future updates to the binary wrap binding file.
*/
const float kWrapFileVersion = 1.0f;

const char* CVWrapCmd::kName = "cvWrap";
const char* CVWrapCmd::kNameFlagShort = "-n";
const char* CVWrapCmd::kNameFlagLong = "-name";
const char* CVWrapCmd::kRadiusFlagShort = "-r";
const char* CVWrapCmd::kRadiusFlagLong = "-radius";
const char* CVWrapCmd::kNewBindMeshFlagShort = "-nbm";
const char* CVWrapCmd::kNewBindMeshFlagLong = "-newBindMesh";
const char* CVWrapCmd::kExportFlagShort = "-ex";
const char* CVWrapCmd::kExportFlagLong = "-export";
const char* CVWrapCmd::kImportFlagShort = "-im";
const char* CVWrapCmd::kImportFlagLong = "-import";
const char* CVWrapCmd::kBindingFlagShort = "-b";
const char* CVWrapCmd::kBindingFlagLong = "-binding";
const char* CVWrapCmd::kRebindFlagShort = "-rb";
const char* CVWrapCmd::kRebindFlagLong = "-rebind";
const char* CVWrapCmd::kHelpFlagShort = "-h";
const char* CVWrapCmd::kHelpFlagLong = "-help";

/**
  Displays command instructions.
*/
void DisplayHelp() {
  MString help;
  help += "Flags:\n"; 
  help += "-name (-n):          String     Name of the wrap node to create.\n"; 
  help += "-radius (-r):        Double     Sample radius.  Default is 0.1.  The greater the radius,\n"; 
  help += "                                the smoother the deformation but slower performance.\n";
  help += "-newBindMesh (-nbm)  N/A        Creates a new bind mesh, otherwise the existing bind mesh will be used.\n";
  help += "-export (-ex):       String     Path to a file to export the binding to.\n"; 
  help += "-import (-im):       String     Path to a file to import the binding from.\n"; 
  help += "-binding (-b):       String     Path to a file to import the binding from on creation.\n"; 
  help += "-rebind (-rb):       String     The name of the wrap node we are rebinding.\n"; 
  help += "-help (-h)           N/A        Display this text.\n";
  MGlobal::displayInfo(help);
}


CVWrapCmd::CVWrapCmd()
    : radius_(0.1),
      name_("cvWrap#"),
      command_(kCommandCreate),
      useBinding_(false),
      newBindMesh_(false) {
}


MSyntax CVWrapCmd::newSyntax() {
  MSyntax syntax;
  syntax.addFlag(kNameFlagShort, kNameFlagLong, MSyntax::kString);
  syntax.addFlag(kRadiusFlagShort, kRadiusFlagLong, MSyntax::kDouble);
  syntax.addFlag(kNewBindMeshFlagShort, kNewBindMeshFlagLong);
  syntax.addFlag(kExportFlagShort, kExportFlagLong, MSyntax::kString);
  syntax.addFlag(kImportFlagShort, kImportFlagLong, MSyntax::kString);
  syntax.addFlag(kBindingFlagShort, kBindingFlagLong, MSyntax::kString);
  syntax.addFlag(kRebindFlagShort, kRebindFlagLong, MSyntax::kString);
  syntax.addFlag(kHelpFlagShort, kHelpFlagLong);
  syntax.setObjectType(MSyntax::kSelectionList, 0, 255);
  syntax.useSelectionAsDefault(true);
  return syntax;
}


void* CVWrapCmd::creator() {                                
  return new CVWrapCmd;                    
}    


bool CVWrapCmd::isUndoable() const {
  return command_ == kCommandCreate;  // Only creation will be undoable
}


MStatus CVWrapCmd::doIt(const MArgList& args) {
  MStatus status;
    
  status = GatherCommandArguments(args);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  if (command_ == kCommandImport || command_ == kCommandExport) {
    // In import/export mode, get the selected wrap deformer node so we can read/write
    // data from it.
    status = selectionList_.getDependNode(0, oWrapNode_);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MFnDependencyNode fnNode(oWrapNode_);
    if (fnNode.typeId() != CVWrap::id) {
      MGlobal::displayError("No wrap node specified.");
      return MS::kFailure;
    }
  } else if (command_ == kCommandRebind) {
    status = GetGeometryPaths();
    CHECK_MSTATUS_AND_RETURN_IT(status);
    status = Rebind();
    CHECK_MSTATUS_AND_RETURN_IT(status);
  } else {
    // Otherwise get the driver and driven geometry paths.
    status = GetGeometryPaths();
    CHECK_MSTATUS_AND_RETURN_IT(status);

    // Add the cvWrap creation command to the modifier.
    MString command = "deformer -type cvWrap -n \"" + name_ + "\"";
    for (unsigned int i = 0; i < pathDriven_.length(); ++i) {
      MFnDagNode fnDriven(pathDriven_[i]);
      command += " " + fnDriven.partialPathName();
    }
    status = dgMod_.commandToExecute(command);
    CHECK_MSTATUS_AND_RETURN_IT(status);
  }

  return redoIt();
}


MStatus CVWrapCmd::GatherCommandArguments(const MArgList& args) {
  MStatus status;
  MArgDatabase argData(syntax(), args);
  argData.getObjects(selectionList_);
  if (argData.isFlagSet(kHelpFlagShort)) {
    command_ = kCommandHelp;
    DisplayHelp();
    return MS::kSuccess;
  } else if (argData.isFlagSet(kExportFlagShort)) {
    command_ = kCommandExport;
    filePath_ = argData.flagArgumentString(kExportFlagShort, 0, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
  } else if (argData.isFlagSet(kImportFlagShort)) {
    command_ = kCommandImport;
    filePath_ = argData.flagArgumentString(kImportFlagShort, 0, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
  }
  newBindMesh_ = argData.isFlagSet(kNewBindMeshFlagShort);
  if (argData.isFlagSet(kRadiusFlagShort)) {
    radius_ = argData.flagArgumentDouble(kRadiusFlagShort, 0, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    // Make sure radius is positive
    if (radius_ <= 0.0) {
      radius_ = 0.001;
    }
  }
  if (argData.isFlagSet(kNameFlagShort)) {
    name_ = argData.flagArgumentString(kNameFlagShort, 0, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
  }
  if (argData.isFlagSet(kBindingFlagShort)) {
    useBinding_ = true;
    filePath_ = argData.flagArgumentString(kBindingFlagShort, 0, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
  }
  if (argData.isFlagSet(kRebindFlagShort)) {
    command_ = kCommandRebind;
    // Get the specified wrap node to rebind.
    MString wrapNode = argData.flagArgumentString(kRebindFlagShort, 0, &status);
    MSelectionList slist;
    status = slist.add(wrapNode);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    status = slist.getDependNode(0, oWrapNode_);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MFnDependencyNode fnNode(oWrapNode_, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    if (fnNode.typeId() != CVWrap::id) {
      MGlobal::displayError(fnNode.name() + " is not a cvWrap node.");
      return MS::kFailure;
    }
  }
  return MS::kSuccess;
}


MStatus CVWrapCmd::GetGeometryPaths() {
  MStatus status;
  // The driver is selected last
  status = selectionList_.getDagPath(selectionList_.length() - 1, pathDriver_, driverComponents_);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = GetShapeNode(pathDriver_);
  // The driver must be a mesh for this specific algorithm.
  if (!pathDriver_.hasFn(MFn::kMesh)) {
    MGlobal::displayError("cvWrap driver must be a mesh.");
    return MS::kFailure;
  }

  MItSelectionList iter(selectionList_);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  pathDriven_.clear();
  drivenComponents_.clear();
  for (unsigned int i = 0; i < selectionList_.length() - 1; ++i, iter.next()) {
    MDagPath path;
    MObject component;
    iter.getDagPath(path, component);
    status = GetShapeNode(path);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    pathDriven_.append(path);
    drivenComponents_.append(component);
  }
  
  return MS::kSuccess;
}


MStatus CVWrapCmd::redoIt() {
  MStatus status;
  if (command_ == kCommandImport) {
    std::ifstream in(filePath_.asChar(), ios::binary);
    if (!in.is_open()) {
      MGlobal::displayInfo("Unable to open file for importing.");
      CHECK_MSTATUS_AND_RETURN_IT(MS::kFailure);
    }
    BindingIO exporter;
    status = exporter.ImportBinding(in, oWrapNode_);
    in.close();
    CHECK_MSTATUS_AND_RETURN_IT(status);
    return MS::kSuccess;
  } else if (command_ == kCommandExport) {
    std::ofstream out(filePath_.asChar(), ios::binary);
    if (!out.is_open()) {
      MGlobal::displayError("Unable to open file for writing.");
      return MS::kFailure;
    }
    BindingIO exporter;
    status = exporter.ExportBinding(out, oWrapNode_);
    out.close();
    CHECK_MSTATUS_AND_RETURN_IT(status);
    return MS::kSuccess;
  } else if (command_ == kCommandRebind) {
    status = dgMod_.doIt();
    CHECK_MSTATUS_AND_RETURN_IT(status);
    return MS::kSuccess;
  } else if (command_ == kCommandCreate) {
    status = CreateWrapDeformer();
    CHECK_MSTATUS_AND_RETURN_IT(status);
    return MS::kSuccess;
  }
  return MS::kFailure;
}

   
MStatus CVWrapCmd::CreateWrapDeformer() {
  MStatus status;
  // Create the deformer
  status = dgMod_.doIt();
  CHECK_MSTATUS_AND_RETURN_IT(status);
  // Reacquire the paths because on referenced geo, a new driven path is created (the ShapeDeformed).
  status = GetGeometryPaths();
  CHECK_MSTATUS_AND_RETURN_IT(status);
  // Get the created wrap deformer node.
  status = GetLatestWrapNode();
  CHECK_MSTATUS_AND_RETURN_IT(status);

  MFnDependencyNode fnNode(oWrapNode_, &status);
  setResult(fnNode.name());
  CHECK_MSTATUS_AND_RETURN_IT(status);

  // Create a bind mesh so we can run rebind commands.  We need a mesh at the state of the 
  // initial binding in order to properly calculate rebinding information.  We can't use
  // the intermediate mesh for rebinding because we may not be binding at the rest pose.
  // Check if this driver already has a bind mesh.
  MDagPath pathBindMesh;
  status = GetExistingBindMesh(pathBindMesh);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  if (newBindMesh_ || !pathBindMesh.isValid()) {
    // No bind mesh exists or the user wants to force create a new one.
    status = CreateBindMesh(pathBindMesh);
    CHECK_MSTATUS_AND_RETURN_IT(status);
  }
  status = ConnectBindMesh(pathBindMesh);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  if (useBinding_) {
    // Import a pre-existing binding.
    std::ifstream in(filePath_.asChar(), ios::binary);
    if (!in.is_open()) {
      MGlobal::displayInfo("Unable to open file for importing.");
      CHECK_MSTATUS_AND_RETURN_IT(MS::kFailure);
    }
    BindingIO exporter;
    status = exporter.ImportBinding(in, oWrapNode_);
    in.close();
    CHECK_MSTATUS_AND_RETURN_IT(status);
  } else {
    MDGModifier dgMod;
    BindData bindData;
    status = CalculateBinding(pathBindMesh, bindData, dgMod);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    status = dgMod.doIt();
    CHECK_MSTATUS_AND_RETURN_IT(status);
  }

  // Connect the driver mesh to the wrap deformer.
  MFnDagNode fnDriver(pathDriver_);
  MPlug plugDriverMesh = fnDriver.findPlug("worldMesh", false, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = plugDriverMesh.selectAncestorLogicalIndex(0, plugDriverMesh.attribute());
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MPlug plugDriverGeo(oWrapNode_, CVWrap::aDriverGeo);
  MDGModifier dgMod;
  dgMod.connect(plugDriverMesh, plugDriverGeo);
  status = dgMod.doIt();
  CHECK_MSTATUS_AND_RETURN_IT(status);

  return MS::kSuccess;
}


MStatus CVWrapCmd::GetLatestWrapNode() {
  MStatus status;
  MObject oDriven = pathDriven_[0].node();
  
  // Since we use MDGModifier to execute the deformer command, we can't get
  // the created deformer node, so we need to find it in the deformation chain.
  MItDependencyGraph itDG(oDriven,
                          MFn::kGeometryFilt,
                          MItDependencyGraph::kUpstream, 
                          MItDependencyGraph::kDepthFirst,
                          MItDependencyGraph::kNodeLevel, 
                          &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MObject oDeformerNode;
  for (; !itDG.isDone(); itDG.next()) {
    oDeformerNode = itDG.currentItem();
    MFnDependencyNode fnNode(oDeformerNode, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    if (fnNode.typeId() == CVWrap::id) {
      oWrapNode_ = oDeformerNode;
      return MS::kSuccess;
    }
  }
  return MS::kFailure;
}


MStatus CVWrapCmd::CreateBindMesh(MDagPath& pathBindMesh) {
  MStatus status;
  MStringArray duplicate;
  MFnDependencyNode fnWrap(oWrapNode_, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MFnDagNode fnDriver(pathDriver_);

  // Calling mesh.duplicate() can give incorrect results due to tweaks and such.
  // We are doing the duplicate here rather than the MDGModifier because we need the name
  // of the duplicated geometry and it would not be reliable to do it from the modifier.
  MGlobal::executeCommand("duplicate -rr -n " + fnWrap.name() + "Base " + fnDriver.partialPathName(), duplicate);
  status = GetDagPath(duplicate[0], pathBindMesh);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = DeleteIntermediateObjects(pathBindMesh);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  bindMeshes_.append(duplicate[0]);

  // Hide the duplicate
  MFnDagNode fnBindMesh(pathBindMesh, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MPlug plug = fnBindMesh.findPlug("visibility", &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = plug.setBool(false);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  
  return MS::kSuccess;
}


MStatus CVWrapCmd::ConnectBindMesh(MDagPath& pathBindMesh) {
  MStatus status;
  // Connect the bind mesh to the wrap node
  status = GetShapeNode(pathBindMesh);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MFnDagNode fnBindMeshShape(pathBindMesh, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MPlug plugBindMessage = fnBindMeshShape.findPlug("message", false, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MPlug plugBindMesh(oWrapNode_, CVWrap::aBindDriverGeo);
  MDGModifier dgMod;
  dgMod.connect(plugBindMessage, plugBindMesh);
  status = dgMod.doIt();
  CHECK_MSTATUS_AND_RETURN_IT(status);
  return MS::kSuccess;
}


MStatus CVWrapCmd::CalculateBinding(MDagPath& pathBindMesh, BindData& bindData,
                                    MDGModifier& dgMod) {
  MStatus status;
  bindData.radius = radius_;

  // Store the bind mesh information.
  // Pre-gather the data from Maya so we can multithread the binding process
  bindData.driverMatrix = pathBindMesh.inclusiveMatrix();
  MObject oBindMesh = pathBindMesh.node();
  status = bindData.intersector.create(oBindMesh, bindData.driverMatrix);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  // We need the adjacency of each vertex in order to crawl the mesh.
  status = GetAdjacency(pathBindMesh, bindData.adjacency);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MFnMesh fnBindMesh(pathBindMesh, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  fnBindMesh.getPoints(bindData.driverPoints, MSpace::kWorld);
  fnBindMesh.getVertexNormals(false, bindData.driverNormals, MSpace::kWorld);
  bindData.perFaceVertices.resize(fnBindMesh.numPolygons());
  bindData.perFaceTriangleVertices.resize(fnBindMesh.numPolygons());
  MIntArray vertexCount, vertexList, triangleCounts, triangleVertices;
  status = fnBindMesh.getVertices(vertexCount, vertexList);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = fnBindMesh.getTriangles(triangleCounts, triangleVertices);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  for (unsigned int faceId = 0, iter = 0, triIter = 0; faceId < vertexCount.length(); ++faceId) {
    bindData.perFaceVertices[faceId].clear();
    for (int i = 0; i < vertexCount[faceId]; ++i, ++iter) {
      bindData.perFaceVertices[faceId].append(vertexList[iter]);
    }
    bindData.perFaceTriangleVertices[faceId].resize(triangleCounts[faceId]);
    for (int triId = 0; triId < triangleCounts[faceId]; ++triId) {
      bindData.perFaceTriangleVertices[faceId][triId].setLength(3);
      bindData.perFaceTriangleVertices[faceId][triId][0] = triangleVertices[triIter++];
      bindData.perFaceTriangleVertices[faceId][triId][1] = triangleVertices[triIter++];
      bindData.perFaceTriangleVertices[faceId][triId][2] = triangleVertices[triIter++];
    }
  }

  // Calculate the binding for each deformed geometry
  MPlug plugBindData(oWrapNode_, CVWrap::aBindData);
  MFnMatrixData fnMatrixData;
  for (unsigned int geomIndex = 0; geomIndex < pathDriven_.length(); ++geomIndex) {
    // Get the plugs to the binding attributes for this geometry
    MPlug plugBind = plugBindData.elementByLogicalIndex(geomIndex, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MPlug plugSampleWeights = plugBind.child(CVWrap::aSampleWeights, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MPlug plugSampleVerts = plugBind.child(CVWrap::aSampleComponents, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MPlug plugSampleBindMatrix = plugBind.child(CVWrap::aBindMatrix, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MPlug plugTriangleVerts = plugBind.child(CVWrap::aTriangleVerts, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MPlug plugBarycentricWeights = plugBind.child(CVWrap::aBarycentricWeights, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);

    // Use the intermediate object for the binding.  This assumes the intermediate object
    // has the same component count as the displayed shape.
    MDagPath pathDriven(pathDriven_[geomIndex]);
    status = GetShapeNode(pathDriven, true);
    if (MFAIL(status)) {
      pathDriven = pathDriven_[geomIndex];
    }
    MItGeometry itGeo(pathDriven, drivenComponents_[geomIndex], &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    int geoCount = itGeo.count();

    status = itGeo.allPositions(bindData.inputPoints, MSpace::kWorld);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    bindData.sampleIds.resize(itGeo.count());
    bindData.weights.resize(itGeo.count());
    bindData.bindMatrices.setLength(itGeo.count());
    bindData.coords.resize(itGeo.count());
    bindData.triangleVertices.resize(itGeo.count());

    // Send off the threads to calculate the binding.
    ThreadData<BindData> threadData[TASK_COUNT];
    CreateThreadData<BindData>(TASK_COUNT, itGeo.count(), &bindData, threadData);
    MThreadPool::init();
    MThreadPool::newParallelRegion(CreateTasks, (void *)threadData);
    MThreadPool::release();

    for (int ii = 0; !itGeo.isDone(); itGeo.next(), ++ii) {
      // Store all the binding data for this component
      // Note for nurbs surfaces the indices may not be continuous.
      int logicalIndex = itGeo.index();
      // Store sample vert ids.
      MFnIntArrayData fnIntData;
      MObject oIntData = fnIntData.create(bindData.sampleIds[ii], &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      MPlug plugSampleVertsElement = plugSampleVerts.elementByLogicalIndex(logicalIndex, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      status = dgMod.newPlugValue(plugSampleVertsElement, oIntData);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      // Store sample weights
      MFnDoubleArrayData fnDoubleData;
      MObject oDoubleData = fnDoubleData.create(bindData.weights[ii], &status);
      assert(bindData.weights[ii].length() > 0);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      MPlug plugSampleWeightsElement = plugSampleWeights.elementByLogicalIndex(logicalIndex,
                                                                               &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      status = dgMod.newPlugValue(plugSampleWeightsElement, oDoubleData);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      // Store bind matrix
      MObject oMatrixData = fnMatrixData.create(bindData.bindMatrices[ii], &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      MPlug plugSampleBindMatrixElement = plugSampleBindMatrix.elementByLogicalIndex(logicalIndex,
                                                                                     &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      status = dgMod.newPlugValue(plugSampleBindMatrixElement, oMatrixData);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      // Store triangle vertices
      MFnNumericData fnNumericData;
      MObject oNumericData = fnNumericData.create(MFnNumericData::k3Int, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      status = fnNumericData.setData3Int(bindData.triangleVertices[ii][0],
                                         bindData.triangleVertices[ii][1],
                                         bindData.triangleVertices[ii][2]);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      MPlug plugTriangleVertsElement = plugTriangleVerts.elementByLogicalIndex(logicalIndex,
                                                                               &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      status = dgMod.newPlugValue(plugTriangleVertsElement, oNumericData);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      // Store barycentric coordinates
      oNumericData = fnNumericData.create(MFnNumericData::k3Float, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      status = fnNumericData.setData3Float(bindData.coords[ii][0], bindData.coords[ii][1],
                                           bindData.coords[ii][2]);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      MPlug plugBarycentricWeightsElement = plugBarycentricWeights.elementByLogicalIndex(
        logicalIndex, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      status = dgMod.newPlugValue(plugBarycentricWeightsElement, oNumericData);
      CHECK_MSTATUS_AND_RETURN_IT(status);
    }
  }
  return MS::kSuccess;
}


void CVWrapCmd::CreateTasks(void *data, MThreadRootTask *pRoot) {
  ThreadData<BindData>* threadData = static_cast<ThreadData<BindData>*>(data);

  if (threadData) {
    int numTasks = threadData[0].numTasks;
    for(int i = 0; i < numTasks; i++) {
      MThreadPool::createTask(CalculateBindingTask, (void *)&threadData[i], pRoot);
    }
    MThreadPool::executeAndJoin(pRoot);
  }
}

bool SortCoords(std::pair<int, float> lhs, std::pair<int, float> rhs) {
  return (lhs.second > rhs.second); 
}


MThreadRetVal CVWrapCmd::CalculateBindingTask(void *pParam) {
  ThreadData<BindData>* pThreadData = static_cast<ThreadData<BindData>*>(pParam);
  double*& alignedStorage = pThreadData->alignedStorage;
  BindData* pData = pThreadData->pData;
  MMeshIntersector& intersector = pData->intersector;
  MMeshIntersector& subsetIntersector = pData->subsetIntersector;
  MPointArray& inputPoints = pData->inputPoints;
  MPointArray& driverPoints = pData->driverPoints;
  MFloatVectorArray& driverNormals = pData->driverNormals;
  std::vector<std::set<int> >& adjacency = pData->adjacency;
  std::vector<MIntArray>& sampleIds = pData->sampleIds;
  std::vector<MDoubleArray>& weights = pData->weights;
  std::vector<BaryCoords>& coords = pData->coords;
  std::vector<MIntArray>& triangleVertices = pData->triangleVertices;
  MMatrixArray& bindMatrices = pData->bindMatrices;

  double radius = pData->radius;

  MMatrix& driverMatrix = pData->driverMatrix;
  std::vector<MIntArray>& perFaceVertices = pData->perFaceVertices;
  std::vector<std::vector<MIntArray> >& perFaceTriangleVertices  = pData->perFaceTriangleVertices;

  unsigned int taskStart = pThreadData->start;
  unsigned int taskEnd = pThreadData->end;

  // Pre-allocate the aligned storage for intrinsics calculation so we are not dynamically allocating
  // memory in the loop.
  std::vector<std::pair<int, float> > sortedCoords(3);
  for (unsigned int i = taskStart; i < taskEnd; ++i) {
    if (i >= inputPoints.length()) {
      break;
    }
    // We need to calculate a bind matrix for each component.
    // The closest point will be the origin of the coordinate system.
    // The weighted normal of the vertices in the sample radius will be one axis.
    // The weight vector from the closest point to the sample vertices will be the other axis.

    MPoint inputPoint = inputPoints[i];
    MPointOnMesh pointOnMesh;
    if (subsetIntersector.isCreated()) {
      // If we are rebinding, limit the closest point to the subset.
      subsetIntersector.getClosestPoint(inputPoint, pointOnMesh);
      inputPoint = MPoint(pointOnMesh.getPoint()) * driverMatrix;
    }

    intersector.getClosestPoint(inputPoint, pointOnMesh);
    int faceId = pointOnMesh.faceIndex();
    int triangleId = pointOnMesh.triangleIndex();

    // Put point in world space so we can calculate the proper bind matrix.
    MPoint closestPoint = MPoint(pointOnMesh.getPoint()) * driverMatrix;

    // Get barycentric coordinates of closestPoint
    triangleVertices[i] = perFaceTriangleVertices[faceId][triangleId];
    GetBarycentricCoordinates(closestPoint, driverPoints[triangleVertices[i][0]],
                              driverPoints[triangleVertices[i][1]],
                              driverPoints[triangleVertices[i][2]],
                              coords[i]);

    // Sort coords highest to lowest so we can easility calculate the up vector
    for (int j = 0; j < 3; ++j) {
      sortedCoords[j] = std::pair<int, float>(triangleVertices[i][j], coords[i][j]);
    }
    std::sort(sortedCoords.begin(), sortedCoords.end(), SortCoords);
    for (int j = 0; j < 3; ++j) {
      triangleVertices[i][j] = sortedCoords[j].first;
      coords[i][j] = sortedCoords[j].second;
    }

    // Get vertices of closest face so we can crawl out from them.
    MIntArray& vertexList = perFaceVertices[faceId];

    // Crawl the surface to find all the vertices within the sample radius.
    std::map<int, double> distances;
    CrawlSurface(closestPoint, vertexList, driverPoints, radius, adjacency, distances);

    // Calculate the weight values per sampled vertex
    CalculateSampleWeights(distances, radius, sampleIds[i], weights[i]);

    // Get the components that form the orthonormal basis.
    MPoint origin;
    MVector up;
    MVector normal;
    CalculateBasisComponents(weights[i], coords[i], triangleVertices[i], driverPoints,
                             driverNormals, sampleIds[i], alignedStorage, origin, up, normal);
    CreateMatrix(origin, normal, up, bindMatrices[i]);
    bindMatrices[i] = bindMatrices[i].inverse();
  }
  return 0;
}


MStatus CVWrapCmd::GetExistingBindMesh(MDagPath &pathBindMesh) {
  MStatus status;
  MObject oDriver = pathDriver_.node();
  MFnDependencyNode fnDriver(oDriver, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  
  // We'll find the bind mesh associated with the driver mesh by traversing the mesh connections
  // through the cvWrap node.
  MPlug plugOutGeom = fnDriver.findPlug("worldMesh", false, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = plugOutGeom.selectAncestorLogicalIndex(0, plugOutGeom.attribute());
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MPlugArray geomPlugs;
  plugOutGeom.connectedTo(geomPlugs, false, true);
  for (unsigned int i = 0; i < geomPlugs.length(); i++) {
    // First iterate through the outMesh connections to find a cvWrap node.
    MObject oThisNode = geomPlugs[i].node();
    MFnDependencyNode fnNode(oThisNode);
    if (fnNode.typeId() == CVWrap::id) {
      status = GetBindMesh(oThisNode, pathBindMesh);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      return MS::kSuccess;
    }
  }
  return MS::kSuccess;
}


MStatus CVWrapCmd::Rebind() {
  MStatus status;

  // Create bind mesh based off of specified faces
  MDagPath pathDriverSubset;
  status = CreateRebindSubsetMesh(pathDriverSubset);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  // Initialize the subset intersector to enable the rebind during the threaded calculation.
  BindData bindData;
  MObject oBindSubsetMesh = pathDriverSubset.node();
  status = bindData.subsetIntersector.create(oBindSubsetMesh, pathDriverSubset.inclusiveMatrix());
  CHECK_MSTATUS_AND_RETURN_IT(status);

  MDagPath pathBindMesh;
  status = GetBindMesh(oWrapNode_, pathBindMesh);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  status = CalculateBinding(pathBindMesh, bindData, dgMod_);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  // Delete the subset mesh since we don't need it anymore
  pathDriverSubset.pop();
  status = MGlobal::executeCommand("delete " + pathDriverSubset.partialPathName());
  CHECK_MSTATUS_AND_RETURN_IT(status);

  return MS::kSuccess;
}


MStatus CVWrapCmd::GetBindMesh(MObject& oWrapNode, MDagPath& pathBindMesh) {
  MStatus status;
  // Get the bind mesh connected to the message attribute of the wrap deformer
  MPlug plugBindMesh(oWrapNode, CVWrap::aBindDriverGeo);
  MPlugArray plugs;
  plugBindMesh.connectedTo(plugs, true, false, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  if (plugs.length() == 0) {
    MGlobal::displayError("Unable to rebind.  No bind mesh is connected.");
    return MS::kFailure;
  }
  MObject oBindMesh = plugs[0].node();
  status = MDagPath::getAPathTo(oBindMesh, pathBindMesh);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  return MS::kSuccess;
}



MStatus CVWrapCmd::CreateRebindSubsetMesh(MDagPath& pathDriverSubset) {
  // We will create the mesh subset by deleting all the non-selected faces.
  MStatus status;

  MDagPath pathBindMesh;
  status = GetBindMesh(oWrapNode_, pathBindMesh);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MFnMesh fnBindMesh(pathBindMesh);

  // Duplicate the bind mesh to create subset
  MStringArray duplicate;
  // Calling mesh.duplicate() gave jacked results.
  status = MGlobal::executeCommand("duplicate -rr " + fnBindMesh.partialPathName(), duplicate);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = GetDagPath(duplicate[0], pathDriverSubset);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = DeleteIntermediateObjects(pathDriverSubset);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  // Get selected driver faces
  MFnSingleIndexedComponent fnDriverComp(driverComponents_, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MIntArray driverFaces;
  status = fnDriverComp.getElements(driverFaces);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  int numFacesToDelete = fnBindMesh.numPolygons() - driverFaces.length();
  if (numFacesToDelete) {
    // Get all the face ids to delete.
    MIntArray facesToDelete;
    int selectedFaceIndex = 0;
    for (int i = 0; i < fnBindMesh.numPolygons(); i++) {
      if (i != driverFaces[selectedFaceIndex]) {
        facesToDelete.append(i);
      } else {
        selectedFaceIndex++;
      }
    }

    MFnSingleIndexedComponent fnDeleteComp;
    MObject oFacesToDelete = fnDeleteComp.create(MFn::kMeshPolygonComponent, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    status = fnDeleteComp.addElements(facesToDelete);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MSelectionList deleteList;
    status = deleteList.add(pathDriverSubset, oFacesToDelete);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    status = MGlobal::setActiveSelectionList(deleteList);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    status = MGlobal::executeCommand("delete;");
    CHECK_MSTATUS_AND_RETURN_IT(status);
    // Reacquire the the dag path since it is invalid now after deleting the faces.
    status = GetDagPath(duplicate[0], pathDriverSubset);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    status = GetShapeNode(pathDriverSubset);
    CHECK_MSTATUS_AND_RETURN_IT(status);
  }
  return MS::kSuccess;
}


MStatus CVWrapCmd::undoIt() {
  MStatus status;
  status = dgMod_.undoIt();
  CHECK_MSTATUS_AND_RETURN_IT(status);

  if (bindMeshes_.length()) {
    // Delete any created bind meshes.
    MDGModifier mod;
    for (unsigned int i = 0; i < bindMeshes_.length(); i++) {
      status = mod.commandToExecute("delete " + bindMeshes_[i]);
      CHECK_MSTATUS_AND_RETURN_IT(status);
    }
    status = mod.doIt();
    CHECK_MSTATUS_AND_RETURN_IT(status);
    bindMeshes_.clear();
  }

  return MS::kSuccess;
}
