#include "cvWrapCmd.h"
#include "cvWrapDeformer.h"
#include "common.h"

#include <maya/MArgDatabase.h>
#include <maya/MFnMatrixData.h>
#include <maya/MItDependencyGraph.h>
#include <maya/MItSelectionList.h>
#include <maya/MMeshIntersector.h>
#include <maya/MSyntax.h>


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
const char* CVWrapCmd::kUVSetFlagShort = "-uv";
const char* CVWrapCmd::kUVSetFlagLong = "-uvSet";
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
  help += "-uvSet (-uv):        String     The UV set to use in the bind process.\n"; 
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
  syntax.addFlag(kHelpFlagShort, kHelpFlagLong);
  syntax.setObjectType(MSyntax::kSelectionList, 0, 255);
  syntax.useSelectionAsDefault(true);
  return syntax;
}


void* CVWrapCmd::creator() {                                
  return new CVWrapCmd;                    
}    


bool CVWrapCmd::isUndoable() const {
  return true;
}


MStatus CVWrapCmd::doIt(const MArgList& args) {
  MStatus status;
    
  GatherCommandArguments(args);

  if (command_ == kCommandImport && command_ == kCommandExport) {
    // In import/export mode, get the selected wrap deformer node so we can read/write
    // data from it.
    status = selectionList_.getDependNode(0, oWrapNode_);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MFnDependencyNode fnNode(oWrapNode_);
    if (fnNode.typeId() != CVWrap::id) {
      MGlobal::displayError("No wrap node specified.");
      return MS::kFailure;
    }
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
  if (argData.isFlagSet(kUVSetFlagShort)) {
    uvset_ = argData.flagArgumentString(kUVSetFlagShort, 0, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
  }
  argData.getObjects(selectionList_);
  return MS::kSuccess;
}


MStatus CVWrapCmd::GetGeometryPaths() {
  MStatus status;
  status = selectionList_.getDagPath(0, pathDriver_);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = GetShapeNode(pathDriver_);
  // The driver must be a mesh for this specific algorithm.
  if (!pathDriver_.hasFn(MFn::kMesh)) {
    MGlobal::displayError("cvWrap driver must be a mesh.");
    return MS::kFailure;
  }
  MItSelectionList iter(selectionList_);
  iter.next(); // Skip the first selected mesh which is the driver
  CHECK_MSTATUS_AND_RETURN_IT(status);
  pathDriven_.clear();
  for ( ; !iter.isDone(); iter.next()) {
      MDagPath path;
      MObject component;
      iter.getDagPath(path, component);
      status = GetShapeNode(path);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      pathDriven_.append(path);
  }
  
  return MS::kSuccess;
}


MStatus CVWrapCmd::redoIt() {
  MStatus status;
  if (command_ == kCommandImport) {
    status = ImportBinding();
    CHECK_MSTATUS_AND_RETURN_IT(status);
    return MS::kSuccess;
  } else if (command_ == kCommandExport) {
    status = ExportBinding();
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

  if (useBinding_) {
    // Import a pre-existing binding.
    status = ImportBinding();
    CHECK_MSTATUS_AND_RETURN_IT(status);
  } else {
    status = CalculateBinding();
    CHECK_MSTATUS_AND_RETURN_IT(status);
  }

  // Connect the driver mesh to the wrap deformer.
  MFnDagNode fnDriver(pathDriver_);
  MPlug plugDriverMesh = fnDriver.findPlug("worldMesh", false, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = plugDriverMesh.selectAncestorLogicalIndex(0, plugDriverMesh.attribute());
  CHECK_MSTATUS_AND_RETURN_IT(status);

  MPlug plugDriverGeo(oWrapNode_, CVWrap::aDriverGeo);
  MPlug plugBindMesh(oWrapNode_, CVWrap::aBindDriverGeo);

  // Create a bind mesh so we can run rebind commands.  We need a mesh at the state of the 
  // initial binding in order to properly calculate rebinding information.  We can't use
  // the intermediate mesh for rebinding because we may not be binding at the rest pose.

  // Check if this driver already has a bind mesh.
  MDagPath pathBindMesh;
  status = GetExistingBindMesh(pathBindMesh);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  if (newBindMesh_ || !pathBindMesh.isValid()) {
    // No bind mesh exists or the user wants to force create a new one.
    MStringArray duplicate;
    // Calling mesh.duplicate() can give incorrect results.
    // We are doing the duplicate here rather than the MDGModifier because we need the name
    // of the duplicated geometry and it would not be reliable to do it from the modifier.
    MGlobal::executeCommand("duplicate -rr -n " + fnNode.name() + "Base " + fnDriver.partialPathName(), duplicate);
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
  }
  status = GetShapeNode(pathBindMesh);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MFnDagNode fnBindMesh(pathBindMesh, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MPlug plugBindMessage = fnBindMesh.findPlug("message", false, &status);   
  CHECK_MSTATUS_AND_RETURN_IT(status);

  // Connect the driver and bind meshes to the wrap deformer
  MDGModifier dgMod;
  dgMod.connect(plugDriverMesh, plugDriverGeo);
  dgMod.connect(plugBindMessage, plugBindMesh);
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


MStatus CVWrapCmd::CalculateBinding() {
  MStatus status;
  MMeshIntersector intersector;
  MObject oDriver = pathDriver_.node();
  MMatrix driverMatrix = pathDriver_.inclusiveMatrix();
  status = intersector.create(oDriver, driverMatrix);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  MFnMesh fnDriverMesh(pathDriver_, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  bool uvSetExists = false;
  if (uvset_.length()) {
    // Make sure uv set exists
    MStringArray uvSets;
    status = fnDriverMesh.getUVSetNames(uvSets);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    for (unsigned int i = 0; i < uvSets.length(); ++i) {
      if (uvset_ == uvSets[i]) {
        uvSetExists = true;
        break;
      }
    }
    if (!uvSetExists) {
      // Fall back to the current UV set
      MString message = "UV set " + uvset_ + " does not exist on " + fnDriverMesh.name() +
                        ", defaulting to current uv set.";
      MGlobal::displayWarning(message);
    }
  }
  if (uvset_.length() == 0 || !uvSetExists) {
    // If no UV set was specified use the current UV set.
    status = fnDriverMesh.getCurrentUVSetName(uvset_);
    CHECK_MSTATUS_AND_RETURN_IT(status);
  }

  // We need the adjacency of each vertex in order to crawl the mesh.
  std::vector<std::set<int> > adjacency;
  status = GetAdjacency(pathDriver_, adjacency);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  MPointArray points;
  MFloatVectorArray normals, tangents;
  fnDriverMesh.getPoints(points, MSpace::kWorld);
  fnDriverMesh.getVertexNormals(false, normals, MSpace::kWorld);
  fnDriverMesh.getTangents(tangents, MSpace::kWorld);

  MFnSingleIndexedComponent fnComp;
  MFnComponentListData fnCompData;
  
  MPlug plugBindData(oWrapNode_, CVWrap::aBindData);
  
  MFnNumericData fnNumericData;
  MFnMatrixData fnMatrixData;
  for (unsigned int geomIndex = 0; geomIndex < pathDriven_.length; ++geomIndex) {
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
    MItGeometry itGeo(pathDriven_[geomIndex], &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);

    // Start progress bar
    StartProgress("Binding wrap...", itGeo.count());

    for (int ii = 0; !itGeo.isDone(); itGeo.next(), ++ii) {
      // We need to calculate a bind matrix for each component.
      // The closest point will be the origin of the coordinate system.
      // The weighted normal of the vertices in the sample radius will be one axis.
      // The tangent at the closest point will be another axis.

      // Get the closest point and faceId.  The close
      MPoint inputPoint = itGeo.position(MSpace::kWorld, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      MPointOnMesh pointOnMesh;
      status = intersector.getClosestPoint(inputPoint, pointOnMesh);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      int faceId = pointOnMesh.faceIndex();
      int triangleId = pointOnMesh.triangleIndex();
      // Put point in world space so we can calculate the proper bind matrix.
      MPoint closestPoint = MPoint(pointOnMesh.getPoint()) * driverMatrix;

      // Get vertices of closest face so we can crawl out from them.
      MIntArray vertexList;
      status = fnDriverMesh.getPolygonVertices(faceId, vertexList);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      // Crawl the surface to find all the vertices within the sample radius.
      std::map<int, double> distances;
      vertexList.clear();
      CrawlSurface(closestPoint, vertexList, points, radius_, adjacency, distances);

      // Calculate the weight values per sampled vertex
      MIntArray sampleIds;
      MDoubleArray weights;
      CalculateSampleWeights(distances, radius_, sampleIds, weights);

      // Get barycentric coordinates of closestPoint
      int triangleVertices[3];
      status = fnDriverMesh.getPolygonTriangleVertices(faceId, triangleId, triangleVertices);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      BaryCoords coords;
      GetBarycentricCoordinates(closestPoint, points[triangleVertices[0]],
                                points[triangleVertices[1]], points[triangleVertices[2]],
                                coords);

      // Calculate the up vector from the tangents
      MVector up;
      for (int j = 0; j < 3; ++j) {
        int tangentId = fnDriverMesh.getTangentId(faceId, triangleVertices[j]);
        up += tangents[tangentId] * coords[j];
      }
      up.normalize();

      // Calculate the weighted normal of the coordinate system
      MVector normal;
      for (unsigned int j = 0; j < weights.length(); j++) {
        normal += MVector(normals[sampleIds[j]]) * weights[j];
      }
      normal.normalize();

      // Store all the binding data for this component
      // Note for nurbs surfaces the indices may not be continuous.
      int logicalIndex = itGeo.index();

      // Store sample vert ids.
      MObject oComp = fnComp.create(MFn::kMeshVertComponent, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      status = fnComp.addElements(sampleIds);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      MObject oCompData = fnCompData.create();
      fnCompData.add(oComp);
      plugSampleVerts.elementByLogicalIndex(logicalIndex, &status).setMObject(oCompData);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      // Store sample weights
      MFnDoubleArrayData fnDoubleData;
      MObject oDoubleData = fnDoubleData.create(weights, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      plugSampleWeights.elementByLogicalIndex(logicalIndex, &status).setMObject(oDoubleData);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      // Store bind matrix
      MMatrix matrix;
      CreateMatrix(closestPoint, normal, up, matrix);
      matrix = matrix.inverse();
      MObject oMatrixData = fnMatrixData.create(matrix, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      plugSampleBindMatrix.elementByLogicalIndex(logicalIndex, &status).setMObject(oMatrixData);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      // Store triangle vertices
      MObject oNumericData = fnNumericData.create(MFnNumericData::k3Int, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      status = fnNumericData.setData3Int(triangleVertices[0], triangleVertices[1], triangleVertices[2]);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      plugTriangleVerts.elementByLogicalIndex(logicalIndex, &status).setMObject(oNumericData);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      // Store barycentric coordinates
      oNumericData = fnNumericData.create(MFnNumericData::k3Float, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      status = fnNumericData.setData3Float(coords[0], coords[1], coords[2]);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      plugBarycentricWeights.elementByLogicalIndex(logicalIndex, &status).setMObject(oNumericData);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      // Increment the progress bar on every 250 verts because updating the UI is slow.
      if (ii % 250 == 0 && ii != 0) {
        StepProgress(250);
        if (ProgressCancelled()) {
          break;
        }
      }
    }

    EndProgress();
  }
  return MS::kSuccess;
}


MStatus CVWrapCmd::GetExistingBindMesh(MDagPath &pathBindMesh) {
  MStatus status;
  MObject oDriver = pathDriver_.node();
  MFnDependencyNode fnDriver(oDriver, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  
  // We'll find the bind mesh associated with the driver mesh by traversing the mesh connections
  // through the cvWrap node.
  MPlug plugOutGeom = fnDriver.findPlug("outMesh", false, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MPlugArray geomPlugs;
  plugOutGeom.connectedTo(geomPlugs, false, true);
  for (unsigned int i = 0; i < geomPlugs.length(); i++) {
    // First iterate through the outMesh connections to find a cvWrap node.
    MObject oThisNode = geomPlugs[i].node();
    MFnDependencyNode fnNode(oThisNode);
    if (fnNode.typeId() == CVWrap::id) {
      // Get bind wrap mesh from wrap node
      MPlug plugBindWrapMesh = fnNode.findPlug(CVWrap::aBindDriverGeo, false, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      MPlugArray bindPlugs;
      plugBindWrapMesh.connectedTo(bindPlugs, true, false);
      if (bindPlugs.length() > 0) {
        // If a bind mesh is connected, use it!
        MObject oBindMesh = bindPlugs[0].node();
        MFnDagNode fnBindDag(oBindMesh, &status);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        status = fnBindDag.getPath(pathBindMesh);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        return MS::kSuccess;
      }
    }
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


MStatus CVWrapCmd::ExportBinding()
{
    MStatus status;
    MFnDependencyNode fnWrapNode(oWrapNode_, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);

    if (fnWrapNode.typeId() != cvWrap::id)
    {
        CHECK_MSTATUS_AND_RETURN_IT(MS::kFailure);
    }

    std::ofstream out(filePath_.asChar(), ios::binary);
    if (!out.is_open())
    {
        MGlobal::displayError("Unable to open file for writing.");
        CHECK_MSTATUS_AND_RETURN_IT(MS::kFailure);
    }

    MFnSingleIndexedComponent fnComp;
    MFnComponentListData fnCompData;
    MFnNumericData fnNumericData;
    MFnMatrixData fnMatrixData;
    MMatrix matrix;
    MObject oTemp;
    MFloatArray tempCoords(3);
    MIntArray triangleVerts(3);
    MDoubleArray tempWeights;
    MIntArray sampleIds;
    MIntArray tempIds;
    MObject oComp;

    MPlug plugSampleWeights(oWrapNode_, cvWrap::aSampleWeights);
    MPlug plugSampleVerts(oWrapNode_, cvWrap::aSampleComponents);
    MPlug plugSampleBindMatrix(oWrapNode_, cvWrap::aBindMatrix);
    MPlug plugTriangleVerts(oWrapNode_, cvWrap::aTriangleVerts);
    MPlug plugBarycentricWeights(oWrapNode_, cvWrap::aBarycentricWeights);

    float version = 1.0f;
    out.write((char *)&version, sizeof(float));

    unsigned int numVerts = plugSampleWeights.numElements();
    out.write((char *)(&numVerts), sizeof(int));

    // Start progress bar
    StartProgress("Exporting wrap...", numVerts);
    for (unsigned int i = 0; i < numVerts; ++i) {
        // Export sample vertex ids
        oTemp = plugSampleVerts.elementByPhysicalIndex(i, &status).asMObject();
        CHECK_MSTATUS_AND_RETURN_IT(status);
        status = fnCompData.setObject(oTemp);
        CHECK_MSTATUS_AND_RETURN_IT(status);

        sampleIds.clear();
        for (unsigned int j = 0; j < fnCompData.length(); j++)
        {
            oComp = fnCompData[j];
            if (oComp.hasFn(MFn::kSingleIndexedComponent))
            {
                status = fnComp.setObject(oComp);
                CHECK_MSTATUS_AND_RETURN_IT(status);
                fnComp.getElements(tempIds);
                for(unsigned int k = 0; k < tempIds.length(); k++)
                {
                    sampleIds.append(tempIds[k]);
                }
            }
        }
        writeAttribute(out, sampleIds);

        // Export sample weights
        oTemp = plugSampleWeights.elementByPhysicalIndex(i, &status).asMObject();
        CHECK_MSTATUS_AND_RETURN_IT(status);
        MFnDoubleArrayData fnDoubleData(oTemp);
        tempWeights = fnDoubleData.array();
        writeAttribute(out, tempWeights);

        // Export bind matrix
        oTemp = plugSampleBindMatrix.elementByPhysicalIndex(i, &status).asMObject();
        CHECK_MSTATUS_AND_RETURN_IT(status);
        status = fnMatrixData.setObject(oTemp);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        writeAttribute(out, fnMatrixData.matrix());

        // Export triangle vertices
        oTemp = plugTriangleVerts.elementByPhysicalIndex(i, &status).asMObject();
        CHECK_MSTATUS_AND_RETURN_IT(status);
        status = fnNumericData.setObject(oTemp);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        fnNumericData.getData3Int(triangleVerts[0], triangleVerts[1], triangleVerts[2]);
        writeAttribute(out, triangleVerts);

        // Export barycentric coordinates
        oTemp = plugBarycentricWeights.elementByPhysicalIndex(i, &status).asMObject();
        CHECK_MSTATUS_AND_RETURN_IT(status);
        status = fnNumericData.setObject(oTemp);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        fnNumericData.getData3Float(tempCoords[0], tempCoords[1], tempCoords[2]);
        writeAttribute(out, tempCoords);

        // Increment the progress bar
        if (i % 250 == 0 && i != 0) {
          StepProgress(250);
          if(ProgressCancelled()) {
            break;
          }
        }
    }
    EndProgress();
    out.close();

    MGlobal::displayInfo("Wrap binding exported.");

    return status;
}


MStatus CVWrapCmd::writeAttribute(std::ofstream &out, const MIntArray &attribute)
{
    unsigned int length = attribute.length();
    out.write((char *)(&length), sizeof(int));
    if (length > 0)
    {
        int * pAttr = new int[length];
        attribute.get(pAttr);
        out.write((char *)pAttr, length * sizeof(int));
        delete [] pAttr;
    }
    return MS::kSuccess;
}


MStatus CVWrapCmd::writeAttribute(std::ofstream &out, const MDoubleArray &attribute)
{
    unsigned int length = attribute.length();
    out.write((char *)(&length), sizeof(int));
    if (length > 0)
    {
        double * pAttr = new double[length];
        attribute.get(pAttr);
        out.write((char *)pAttr, length * sizeof(double));
        delete [] pAttr;
    }
    return MS::kSuccess;
}


MStatus CVWrapCmd::writeAttribute(std::ofstream &out, const MFloatArray &attribute)
{
    unsigned int length = attribute.length();
    out.write((char *)(&length), sizeof(int));
    if (length > 0)
    {
        float * pAttr = new float[length];
        attribute.get(pAttr);
        out.write((char *)pAttr, length * sizeof(float));
        delete [] pAttr;
    }
    return MS::kSuccess;
}


MStatus CVWrapCmd::writeAttribute(std::ofstream &out, const MMatrix &attribute)
{
    unsigned int length = 16;
    out.write((char *)(&length), sizeof(int));
    double pAttr[16];
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            pAttr[(i * 4) + j] = attribute[i][j];
        }
    }
    out.write((char *)pAttr, length * sizeof(double));
    return MS::kSuccess;
}


MStatus CVWrapCmd::ImportBinding()
{
    MStatus status;

    std::ifstream in(filePath_.asChar(), ios::binary);
    if (!in.is_open())
    {
        MGlobal::displayInfo("Unable to open file for importing.");
        CHECK_MSTATUS_AND_RETURN_IT(MS::kFailure);
    }

    MFnSingleIndexedComponent fnComp;
    MFnComponentListData fnCompData;
    MFnNumericData fnNumericData;
    MFnMatrixData fnMatrixData;
    MMatrix matrix;
    MObject oTemp, oData;
    MPlug plugTemp;
    MFloatArray tempCoords(3);
    MIntArray triangleVerts(3);
    MDoubleArray tempWeights;
    MIntArray sampleIds;
    MFnDoubleArrayData fnDoubleData;

    MPlug plugSampleWeights(oWrapNode_, cvWrap::aSampleWeights);
    MPlug plugSampleVerts(oWrapNode_, cvWrap::aSampleComponents);
    MPlug plugSampleBindMatrix(oWrapNode_, cvWrap::aBindMatrix);
    MPlug plugTriangleVerts(oWrapNode_, cvWrap::aTriangleVerts);
    MPlug plugBarycentricWeights(oWrapNode_, cvWrap::aBarycentricWeights);

    float version;
    in.read((char *)(&version), sizeof(float));

    int numVerts;
    in.read((char *)(&numVerts), sizeof(int));
    //if (numVerts != numDrivenVerts)
    //{
    //    char buffer[256];
    //    sprintf(buffer, "Mesh has %d vertices. Binding file containts %d vertices", numDrivenVerts, numVerts);
    //    MGlobal::displayError(buffer);
    //    in.close();
    //    return MS::kFailure;
    //}
    StartProgress("Importing wrap...", numVerts);
    for (int i = 0; i < numVerts; i++)
    {
        // Import sample vertex ids
        readAttribute(in, sampleIds);
        oTemp = fnComp.create(MFn::kMeshVertComponent, &status);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        status = fnComp.addElements(sampleIds);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        oData = fnCompData.create(&status);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        status = fnCompData.add(oTemp);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        plugSampleVerts.elementByLogicalIndex(i, &status).setMObject(oData);
        CHECK_MSTATUS_AND_RETURN_IT(status);

        // Import sample weights
        readAttribute(in, tempWeights);
        MFnDoubleArrayData fnDoubleData;
        oData = fnDoubleData.create(tempWeights, &status);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        plugSampleWeights.elementByLogicalIndex(i, &status).setMObject(oData);
        CHECK_MSTATUS_AND_RETURN_IT(status);

        // Import bind matrix
        readAttribute(in, matrix);
        oData = fnMatrixData.create(matrix);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        plugSampleBindMatrix.elementByLogicalIndex(i, &status).setMObject(oData);
        CHECK_MSTATUS_AND_RETURN_IT(status);

        // Import triangle vertices
        readAttribute(in, triangleVerts);
        oData = fnNumericData.create(MFnNumericData::k3Int, &status);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        status = fnNumericData.setData3Int(triangleVerts[0], triangleVerts[1], triangleVerts[2]);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        plugTriangleVerts.elementByLogicalIndex(i, &status).setMObject(oData);
        CHECK_MSTATUS_AND_RETURN_IT(status);

        // Import barycentric coordinates
        readAttribute(in, tempCoords);
        oData = fnNumericData.create(MFnNumericData::k3Float, &status);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        fnNumericData.setData3Float(tempCoords[0], tempCoords[1], tempCoords[2]);
        plugBarycentricWeights.elementByLogicalIndex(i, &status).setMObject(oData);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        // Increment the progress bar
        if (i % 250 == 0 && i != 0) {
          StepProgress(250);
          if (ProgressCancelled()) {
            break;
          }
        }
    }
    EndProgress();
    in.close();
    MGlobal::displayInfo("Wrap binding imported.");

    return MS::kSuccess;
}


MStatus CVWrapCmd::readAttribute(std::ifstream &in, MIntArray &attribute)
{
    attribute.clear();
    int numValues;
    in.read((char *)(&numValues), sizeof(int));
    if (numValues > 0)
    {
        attribute.setLength(numValues);
        int * pValues = new int[numValues];
        in.read((char *)pValues, numValues * sizeof(int));
        for (int i = 0; i < numValues; i++)
        {
            attribute[i] = pValues[i];
        }
        delete [] pValues;
    }

    return MS::kSuccess;
}


MStatus CVWrapCmd::readAttribute(std::ifstream &in, MDoubleArray &attribute)
{
    attribute.clear();
    int numValues;
    in.read((char *)(&numValues), sizeof(int));
    if (numValues > 0)
    {
        attribute.setLength(numValues);
        double * pValues = new double[numValues];
        in.read((char *)pValues, numValues * sizeof(double));
        for (int i = 0; i < numValues; i++)
        {
            attribute[i] = pValues[i];
        }
        delete [] pValues;
    }
    return MS::kSuccess;
}


MStatus CVWrapCmd::readAttribute(std::ifstream &in, MFloatArray &attribute)
{
    attribute.clear();
    int numValues;
    in.read((char *)(&numValues), sizeof(int));
    if (numValues > 0)
    {
        attribute.setLength(numValues);
        float * pValues = new float[numValues];
        in.read((char *)pValues, numValues * sizeof(float));
        for (int i = 0; i < numValues; i++)
        {
            attribute[i] = pValues[i];
        }
        delete [] pValues;
    }
    return MS::kSuccess;
}


MStatus CVWrapCmd::readAttribute(std::ifstream &in, MMatrix &matrix)
{
    int numValues;
    in.read((char *)(&numValues), sizeof(int));
    if (numValues > 0)
    {
        double * pValues = new double[numValues];
        in.read((char *)pValues, numValues * sizeof(double));
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                matrix[i][j] = pValues[(i * 4) + j];
            }
        }
        delete [] pValues;
    }
    return MS::kSuccess;
}



