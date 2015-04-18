#ifndef CVWRAPREBINDCMD_H
#define CVWRAPREBINDCMD_H

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
#include <maya/MMatrixArray.h>
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
#include <maya/MItMeshVertex.h>
#include <maya/MItSelectionList.h>

#include <maya/MFn.h>
#include <maya/MFnAttribute.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnDoubleIndexedComponent.h>
#include <maya/MFnData.h>
#include <maya/MFnMatrixData.h>
#include <maya/MFnMesh.h>
#include <maya/MFnNumericData.h>
#include <maya/MFnNurbsSurface.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MFnSet.h>
#include <maya/MFnSubd.h>
#include <maya/MFnTransform.h>
#include <maya/MFnWeightGeometryFilter.h>


#include <maya/MPxCommand.h>

#include <stdio.h>
#include <math.h>
#include <vector>
#include <map>


class cvWrapRebindCmd : public MPxCommand 
{							
public:															
    cvWrapRebindCmd();							
    ~cvWrapRebindCmd();							
	virtual MStatus	doIt( const MArgList& );					
	virtual MStatus	undoIt();
	virtual MStatus	redoIt();
	virtual bool    isUndoable() const;
	static  void*   creator();		
    static  MSyntax newSyntax();
    
private:
    double          m_radius;
    bool            m_weightMapMultiply;
    MString         m_wrapNode;
    MSelectionList  m_selectionList;
    MIntArray       m_selectedVerts;
    MIntArray       m_selectedFaces;
    MDagPath        m_pathDriver;
    MDagPath        m_pathDriven;
    MObject         m_oCompDriver;
    MObject         m_oCompDriven;
    std::vector<MIntArray>          m_origSampleVerts;
    std::vector<MDoubleArray>       m_origSampleWeights;
    std::vector<MIntArray>          m_origTriangleVertices;
    std::vector<MFloatArray>        m_origBarycentricCoordinates;
    MMatrixArray                    m_origBindMatrices;

    MStatus getPathsToDriverAndDriven(
        MDagPath& pathDriver, MObject& oCompDriver,
        MDagPath& pathDriven, MObject& oCompDriven );
    MStatus getWrapNode( MObject& oWrapNode );
    MStatus getBindMesh( MObject& oWrapNode, MDagPath& pathBindMesh );
    MStatus rebindVertices( MObject& oWrapNode );
    MStatus getConnectedPlug( MPlug& plug, MPlug& connectedPlug );
    MStatus createDriverSubsetMesh( MDagPath& pathDriver, MDagPath& pathSubsetMesh );
    MStatus getDagPath( MString& name, MDagPath& pathNode );
    MStatus getMObject( MString& name, MObject& oNode );
    MStatus deleteIntermediateObjects( MDagPath& pathNode );
    MStatus getIntermediateObject( MDagPath& pathShape, MDagPath& pathIntermediate );

};	

#endif
