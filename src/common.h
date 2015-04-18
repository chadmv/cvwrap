/**
  Contains various helper functions.
*/

#ifndef CVWRAP_COMMON_H
#define CVWRAP_COMMON_H

#include <maya/MDagPath.h>
#include <maya/MPoint.h>
#include <maya/MString.h>

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

double DistanceSquared(const MPoint& p1, const MPoint& p2);

#endif