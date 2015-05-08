#include "cvWrapDeformer.h"
#include "cvWrapCmd.h"

#include <maya/MFnPlugin.h>

MStatus initializePlugin(MObject obj) { 
  MStatus status;
  MFnPlugin plugin(obj, "Chad Vernon", "1.0", "Any");
  status = plugin.registerNode(CVWrap::kName, CVWrap::id, CVWrap::creator, CVWrap::initialize,
                               MPxNode::kDeformerNode);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = plugin.registerCommand(CVWrapCmd::kName, CVWrapCmd::creator, CVWrapCmd::newSyntax);
  CHECK_MSTATUS_AND_RETURN_IT(status);
#if MAYA_API_VERSION >= 201600
	status = MGPUDeformerRegistry::registerGPUDeformerCreator(CVWrap::kName, "cvWrapOverride",
                                                            CVWrapGPU::GetGPUDeformerInfo());
  CHECK_MSTATUS_AND_RETURN_IT(status);
  // Set the load path so we can find the cl kernel.
  CVWrapGPU::pluginLoadPath = plugin.loadPath();
#endif

  return status;
}

MStatus uninitializePlugin( MObject obj) {
  MStatus status;
  MFnPlugin plugin(obj);

#if MAYA_API_VERSION >= 201600
  status = MGPUDeformerRegistry::deregisterGPUDeformerCreator(CVWrap::kName, "cvWrapOverride");
  CHECK_MSTATUS_AND_RETURN_IT(status);
#endif
  status = plugin.deregisterCommand(CVWrapCmd::kName);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = plugin.deregisterNode(CVWrap::id);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  
  return status;
}
