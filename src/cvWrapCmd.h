#ifndef CVWRAPCMD_H
#define CVWRAPCMD_H

#include <maya/MArgList.h>
#include <maya/MDagPath.h>
#include <maya/MDGModifier.h>
#include <maya/MFloatArray.h>
#include <maya/MPlug.h>
#include <maya/MSelectionList.h>
#include <maya/MString.h>
#include <maya/MStringArray.h>

#include <maya/MPxCommand.h>

#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>

/**
  The cvWrap command is used to create new cvWrap deformers and to import and export
  wrap bindings.
*/
class CVWrapCmd : public MPxCommand {							
 public:
  enum CommandMode { kCommandCreate, kCommandExport, kCommandImport, kCommandNewBindMesh,
                     kCommandHelp };
  CVWrapCmd();							
	virtual MStatus	doIt(const MArgList&);
	virtual MStatus	undoIt();
	virtual MStatus	redoIt();
	virtual bool isUndoable() const;
	static void* creator();		
  static MSyntax newSyntax();

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

    static void getBarycentricCoordinates( MPoint& P,
            MPoint& A, MPoint& B, MPoint& C,
            float& a, float& b, float& c );
    static MStatus getAdjacency( MDagPath& pathMesh, std::vector<std::set<int> >& adjacency );
    static MStatus crawlSurface( int vertexIndex,
            MPointArray& points, std::map<int, double>& distances, double sourceDistance,
            MPoint& sourcePoint, double maxDistance, std::vector<std::set<int> >& adjacency );
    static void calculateSampleWeights( std::map<int, double>& distances,
            MIntArray& sampleIds,
            MDoubleArray& weights,
            MDoubleArray& normalizedWeights,
            MPoint& closestPoint,
            MPointArray& points,
            double& totalWeight );

    MStatus getExistingBindMesh( MObject& oDriver, MDagPath &pathBindMesh );
    static MStatus getDagPath( MString& name, MDagPath& pathNode );
    static MStatus deleteIntermediateObjects( MDagPath& pathNode );
    static double distanceSquared( const MPoint& p1, const MPoint& p2 );

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

  MStatus CalculateBinding();


  /**
    Exports the binding information to disk.
  */
  MStatus ExportBinding();
  MStatus writeAttribute( std::ofstream &out, const MIntArray &attribute );
  MStatus writeAttribute( std::ofstream &out, const MDoubleArray &attribute );
  MStatus writeAttribute( std::ofstream &out, const MFloatArray &attribute );
  MStatus writeAttribute( std::ofstream &out, const MMatrix &attribute );

  /**
    Imports the binding information from disk.
  */
  MStatus ImportBinding();
  MStatus readAttribute( std::ifstream &in, MIntArray &attribute );
  MStatus readAttribute( std::ifstream &in, MDoubleArray &attribute );
  MStatus readAttribute( std::ifstream &in, MFloatArray &attribute );
  MStatus readAttribute( std::ifstream &out, MMatrix &attribute );

  MString name_;  /**< Name of cvWrap node to create. */
  double radius_;  /**< Binding sample radius. */
  CommandMode command_;
  MString filePath_;
  bool m_edit;
  bool useBinding_;
  MSelectionList selectionList_;  /**< Selected command input nodes. */
  MObject oWrapNode_;  /**< MObject to the cvWrap node in focus. */
  MDagPath pathDriver_;  /**< Path to the shape wrapping the other shape. */
  MDagPath pathDriven_;  /**< Path to the shape being wrapped. */
  MDGModifier dgMod_;
  MStringArray    m_createdNodes;

  
};	

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

#endif
