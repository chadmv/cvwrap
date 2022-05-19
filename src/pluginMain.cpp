#include "cvWrapDeformer.h"
#include "cvWrapCmd.h"

#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>

PLUGIN_EXPORT MStatus initializePlugin(MObject obj) { 
  MStatus status;
  MFnPlugin plugin(obj, "Chad Vernon", "1.0", "Any");
  status = plugin.registerNode(CVWrap::kName, CVWrap::id, CVWrap::creator, CVWrap::initialize,
                               MPxNode::kDeformerNode);
  CHECK_STATUS_AND_RETURN_IT(status);
  status = plugin.registerCommand(CVWrapCmd::kName, CVWrapCmd::creator, CVWrapCmd::newSyntax);
  CHECK_STATUS_AND_RETURN_IT(status);
#if MAYA_API_VERSION >= 201600
	status = MGPUDeformerRegistry::registerGPUDeformerCreator(CVWrap::kName, "cvWrapOverride",
                                                            CVWrapGPU::GetGPUDeformerInfo());
  CHECK_STATUS_AND_RETURN_IT(status);
  // Set the load path so we can find the cl kernel.
  CVWrapGPU::pluginLoadPath = plugin.loadPath();
#endif


  if (MGlobal::mayaState() == MGlobal::kInteractive) {
    MGlobal::executePythonCommandOnIdle("import cvwrap.menu");
		MGlobal::executePythonCommandOnIdle("cvwrap.menu.create_menuitems()");
  }

  return status;
}

PLUGIN_EXPORT MStatus uninitializePlugin( MObject obj) {
  MStatus status;
  MFnPlugin plugin(obj);

#if MAYA_API_VERSION >= 201600
  status = MGPUDeformerRegistry::deregisterGPUDeformerCreator(CVWrap::kName, "cvWrapOverride");
  CHECK_STATUS_AND_RETURN_IT(status);
#endif
  status = plugin.deregisterCommand(CVWrapCmd::kName);
  CHECK_STATUS_AND_RETURN_IT(status);
  status = plugin.deregisterNode(CVWrap::id);
  CHECK_STATUS_AND_RETURN_IT(status);

  if (MGlobal::mayaState() == MGlobal::kInteractive) {
    MGlobal::executePythonCommandOnIdle("import cvwrap.menu");
		MGlobal::executePythonCommandOnIdle("cvwrap.menu.destroy_menuitems()");
  }
  
  return status;
}
