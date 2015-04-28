#ifndef CVWRAPCMD_H
#define CVWRAPCMD_H

#include <maya/MArgList.h>
#include <maya/MDagPath.h>
#include <maya/MDagPathArray.h>
#include <maya/MDGModifier.h>
#include <maya/MFloatArray.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MMatrixArray.h>
#include <maya/MMeshIntersector.h>
#include <maya/MPlug.h>
#include <maya/MPointArray.h>
#include <maya/MSelectionList.h>
#include <maya/MString.h>
#include <maya/MStringArray.h>
#include <maya/MThreadPool.h>

#include <maya/MPxCommand.h>

#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>

#include "common.h"

struct BindData {
  MPointArray inputPoints;  /**< The world space points of the geometry to be wrapped. */
  MPointArray driverPoints;  /**< The world space points of the driver geometry. */
  MFloatVectorArray driverNormals;  /**< The world space normals of the driver geometry. */
  std::vector<MIntArray> perFaceVertices;  /**< The per-face vertex ids of the driver. */
  std::vector<std::vector<MIntArray> > perFaceTriangleVertices;  /**< The per-face per-triangle vertex ids of the driver. */
  MMeshIntersector intersector;  /**< Closest point intersector on the driver mesh. */
  std::vector<std::set<int> > adjacency;  /**< Driver adjacency for surface crawling. */
  MMatrix driverMatrix;  /**< Driver matrix to convert closest points into world space. */
  double radius;  /**< Max crawl sample radius. */

  /**
    Elements calculated in the threads.
  */
  std::vector<MIntArray> sampleIds;
  std::vector<MDoubleArray> weights;
  MMatrixArray bindMatrices;
  std::vector<BaryCoords> coords;
  std::vector<MIntArray> triangleVertices;
};

/**
  The cvWrap command is used to create new cvWrap deformers and to import and export
  wrap bindings.
*/
class CVWrapCmd : public MPxCommand {              
 public:
  enum CommandMode { kCommandCreate, kCommandExport, kCommandImport, kCommandHelp };
  CVWrapCmd();              
  virtual MStatus  doIt(const MArgList&);
  virtual MStatus  undoIt();
  virtual MStatus  redoIt();
  virtual bool isUndoable() const;
  static void* creator();    
  static MSyntax newSyntax();

  /**
    Distributes the ThreadData objects to the parallel threads.
    @param[in] data The user defined data.  In this case, the ThreadData array.
    @param[in] pRoot Maya's root task.
  */
  static void CreateTasks(void *data, MThreadRootTask *pRoot);
  static MThreadRetVal CalculateBindingTask(void *pParam);

  const static char* kName;  /**< The name of the command. */
  
  /**
    Specifies the name of the cvWrap node.
  */
  const static char* kNameFlagShort;
  const static char* kNameFlagLong;
  
  /**
    Specifies the sample radius of the binding.
  */
  const static char* kRadiusFlagShort;
  const static char* kRadiusFlagLong;

  /**
    Specifies that a new bind mesh should be created.  The bind mesh is only used for rebinding
    vertices and can be deleted at any time.  Sometimes, artists may want to wrap different
    geometry with the same mesh.  By default the command will reuse the same bind mesh for a driver,
    but if new geometry is being wrapped at a different pose, a new bind mesh should be created
    in order to correctly rebind.
  */
  const static char* kNewBindMeshFlagShort;
  const static char* kNewBindMeshFlagLong;

  /**
    Export file path.
  */
  const static char* kExportFlagShort;
  const static char* kExportFlagLong;

  /**
    Import file path.
  */
  const static char* kImportFlagShort;
  const static char* kImportFlagLong;
  
  /**
    Path of a binding on disk rather than calculating binding from scratch.
  */
  const static char* kBindingFlagShort;  
  const static char* kBindingFlagLong;

  /**
    Displays help.
  */
  const static char* kHelpFlagShort;
  const static char* kHelpFlagLong;

 private:
  /**
    Gathers all the command arguments and sets necessary command states.
    @param[in] args Maya MArgList.
  */
  MStatus GatherCommandArguments(const MArgList& args);

  /**
    Acquires the driver and driven dag paths from the input selection list.
  */
  MStatus GetGeometryPaths();

  /**
    Creates a new wrap deformer.
  */
  MStatus CreateWrapDeformer();

  /**
    Gets the latest cvWrap node in the history of the deformed shape.
  */
  MStatus GetLatestWrapNode();

  /**
    Calculates the binding data for the wrap deformer to work.
  */
  MStatus CalculateBinding();
    
  /**
    Gets the MDagPath of any existing bind wrap mesh so we don't have to duplicate it for each
    new wrap.
    @param[out] pathBindMesh Storage for path to an existing bind mesh
  */
  MStatus GetExistingBindMesh(MDagPath &pathBindMesh);

  /**
    Exports the binding information to disk.
  */
  MStatus ExportBinding(std::ofstream& out);

  /**
    Imports the binding information from disk.
  */
  MStatus ImportBinding(std::ifstream& in);

  MString name_;  /**< Name of cvWrap node to create. */
  double radius_;  /**< Binding sample radius. */
  CommandMode command_;
  MString filePath_;
  bool m_edit;
  bool useBinding_;
  bool newBindMesh_;
  MSelectionList selectionList_;  /**< Selected command input nodes. */
  MObject oWrapNode_;  /**< MObject to the cvWrap node in focus. */
  MDagPath pathDriver_;  /**< Path to the shape wrapping the other shape. */
  MDagPathArray pathDriven_;  /**< Paths to the shapes being wrapped. */
  MDGModifier dgMod_;
  MStringArray bindMeshes_;

  
};  


template <typename dataType, typename container>
void WriteAttribute(std::ofstream &out, const container& attribute) {
  unsigned int length = attribute.length();
  out.write((char *)(&length), sizeof(length));
  if (length > 0) {
    dataType * pAttr = new dataType[length];
    attribute.get(pAttr);
    out.write((char *)pAttr, length * sizeof(dataType));
    delete [] pAttr;
  }
}


template <>
void WriteAttribute<double, MMatrix>(std::ofstream &out, const MMatrix& attribute);


template <typename dataType, typename container>
void ReadAttribute(std::ifstream &in, container &attribute) {
  attribute.clear();
  unsigned int length;
  in.read((char *)(&length), sizeof(length));
  if (length > 0) {
    attribute.setLength(length);
    dataType* pValues = new dataType[length];
    in.read((char *)pValues, length * sizeof(int));
    for (unsigned int i = 0; i < length; i++) {
      attribute[i] = pValues[i];
    }
    delete [] pValues;
  }
}

template <>
void ReadAttribute<double, MMatrix>(std::ifstream &in, MMatrix &matrix);


#endif
