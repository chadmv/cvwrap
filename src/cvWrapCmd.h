#ifndef CVWRAPCMD_H
#define CVWRAPCMD_H

#include <maya/MArgDatabase.h>
#include <maya/MArgList.h>
#include <maya/MCommandResult.h>
#include <maya/MDagPath.h>
#include <maya/MDagPathArray.h>
#include <maya/MDGModifier.h>
#include <maya/MDoubleArray.h>
#include <maya/MFloatArray.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MGlobal.h>
#include <maya/MIntArray.h>
#include <maya/MMeshIntersector.h>
#include <maya/MObject.h>
#include <maya/MObjectArray.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MPointArray.h>
#include <maya/MRichSelection.h>
#include <maya/MSelectionList.h>
#include <maya/MString.h>
#include <maya/MStringArray.h>
#include <maya/MSyntax.h>
#include <maya/MString.h>
#include <maya/MDGModifier.h>

#include <maya/MItGeometry.h>
#include <maya/MItDependencyGraph.h>

#include <maya/MPxCommand.h>

#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>


class CVWrapCmd : public MPxCommand {							
 public:															
  CVWrapCmd();							
	virtual MStatus	doIt(const MArgList&);
	virtual MStatus	undoIt();
	virtual MStatus	redoIt();
	virtual bool isUndoable() const;
	static void* creator();		
  static MSyntax newSyntax();

  const static char* kName;

    MStatus getPathsToDriverAndDriven( MDagPath& pathDriver, MDagPath& pathDriven );
    MStatus getLatestWrapNode( MDagPath& pathDriven, MObject& oDeformerNode );
    MStatus calculateBinding( MObject& oWrapNode, MDagPath& pathDriver, MDagPath& pathDriven );
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
    void displayHelp();

private:
    double          m_radius;
    bool            m_newBindMesh;
    bool            m_export;
    bool            m_import;
    bool            m_edit;
    bool            m_displayHelp;
    MString         m_filePath;
    bool            m_useBinding;
    MString         m_name;
    MDagPath        m_pathDriver;
    MDagPath        m_pathDriven;
    MSelectionList  m_selectionList;
    MDGModifier     m_dgMod;
    MStringArray    m_createdNodes;
    MObject         m_oWrapNode;

    MStatus         exportBinding();
    MStatus         writeAttribute( std::ofstream &out, const MIntArray &attribute );
    MStatus         writeAttribute( std::ofstream &out, const MDoubleArray &attribute );
    MStatus         writeAttribute( std::ofstream &out, const MFloatArray &attribute );
    MStatus         writeAttribute( std::ofstream &out, const MMatrix &attribute );

    MStatus         importBinding();
    MStatus         readAttribute( std::ifstream &in, MIntArray &attribute );
    MStatus         readAttribute( std::ifstream &in, MDoubleArray &attribute );
    MStatus         readAttribute( std::ifstream &in, MFloatArray &attribute );
    MStatus         readAttribute( std::ifstream &out, MMatrix &attribute );
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
