#include "common.h"

#include <maya/MGlobal.h>
#include <maya/MFnDagNode.h>
#include <maya/MSelectionList.h>

void StartProgress(const MString& title, unsigned int count) {
  if (MGlobal::mayaState() == MGlobal::kInteractive) {
    MString message = "progressBar -e -bp -ii true -st \"";
    message += title;
    message += "\" -max ";
    message += count;
    message += " $gMainProgressBar;";
    MGlobal::executeCommand(message);
  }
}


void StepProgress(int step) {
  if (MGlobal::mayaState() == MGlobal::kInteractive) {
    MString message = "progressBar -e -s ";
    message += step;
    message += " $gMainProgressBar;";
    MGlobal::executeCommand(message);
  }
}


bool ProgressCancelled() {
  if (MGlobal::mayaState() == MGlobal::kInteractive) {
    int cmdResult = 0;
    MGlobal::executeCommand("progressBar -query -isCancelled $gMainProgressBar", cmdResult);
    return cmdResult != 0;
  }
  return false;
}


void EndProgress() {
  if (MGlobal::mayaState() == MGlobal::kInteractive) {
    MGlobal::executeCommand("progressBar -e -ep $gMainProgressBar;");
  }
}


bool IsShapeNode(MDagPath& path) {
  return path.node().hasFn(MFn::kMesh) ||
         path.node().hasFn(MFn::kNurbsCurve) ||
         path.node().hasFn(MFn::kNurbsSurface);
}


MStatus GetShapeNode(MDagPath& path) {
  MStatus status;

  if (path.hasFn(MFn::kTransform)) {
    unsigned int shapeCount;
    status = path.numberOfShapesDirectlyBelow(shapeCount);
    CHECK_MSTATUS_AND_RETURN_IT(status);

    for (unsigned int i = 0; i < shapeCount; ++i) {
      status = path.extendToShapeDirectlyBelow(i);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      // Make sure it is not an intermediate object.
      MFnDagNode fnNode(path, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      if (!fnNode.isIntermediateObject()) {
        return MS::kSuccess;
      }
      // Go to the next shape if this is an intermediate shape.
      path.pop();
    }
  } else if (IsShapeNode(path)) {
    // This is already a shape node.
    return MS::kSuccess;
  }

  // No valid shape node found.
  return MS::kFailure;
}


MStatus GetDagPath(MString& name, MDagPath& path) {
  MStatus status;
  MSelectionList list;
  status = MGlobal::getSelectionListByName(name, list);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = list.getDagPath(0, path);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  return MS::kSuccess;
}


MStatus DeleteIntermediateObjects(MDagPath& path) {
  MStatus status;
  if (IsShapeNode(path)) {
    path.pop();
  }
  if (path.hasFn(MFn::kTransform)) {
    unsigned int shapeCount;
    status = path.numberOfShapesDirectlyBelow(shapeCount);
    CHECK_MSTATUS_AND_RETURN_IT(status);

    for (unsigned int i = 0; i < shapeCount; ++i) {
      status = path.extendToShapeDirectlyBelow(i);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      MFnDagNode fnNode(path, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      path.pop();
      if (fnNode.isIntermediateObject()) {
        status = MGlobal::executeCommand("delete " + fnNode.partialPathName());
        CHECK_MSTATUS_AND_RETURN_IT(status);
      }
    }
  }
  return MS::kSuccess;
}


double DistanceSquared(const MPoint& p1, const MPoint& p2) {
  double xx = p2.x - p1.x;
  double yy = p2.y - p1.y;
  double zz = p2.z - p1.z;
  return (xx*xx) + (yy*yy) + (zz*zz);
}