#ifndef CVWRAPDEFORMER_H
#define CVWRAPDEFORMER_H

#include <maya/MArrayDataHandle.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MDoubleArray.h>
#include <maya/MFloatArray.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MFnComponentListData.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnData.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnGenericAttribute.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMessageAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnNurbsSurface.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MGlobal.h>
#include <maya/MIntArray.h>
#include <maya/MMatrix.h> 
#include <maya/MMatrixArray.h> 
#include <maya/MNodeMessage.h> 
#include <maya/MPlug.h> 
#include <maya/MPoint.h> 
#include <maya/MPointArray.h> 
#include <maya/MQuaternion.h>
#include <maya/MThreadPool.h>
#include <maya/MTransformationMatrix.h>
#include <maya/MTypeId.h> 
#include <maya/MVector.h>
#include <maya/MVectorArray.h>

#include <maya/MItGeometry.h>
#include <maya/MItMeshVertex.h>

#include <maya/MPxDeformerNode.h>


#include <vector>

struct Sortable
{
    double weight;
    double normalizedWeight;
    int index;
};

typedef struct TaskData 
{
    float envelope;
    MIntArray membership;
    MFloatArray weights;
    MPointArray points;

    MPointArray driverPoints;
    MVectorArray driverNormals;
    std::vector< std::vector <Sortable> > sortedWeights;
    MMatrixArray bindMatrices;
    std::vector<MIntArray> triangleVerts;
    std::vector<MFloatArray> barycentricWeights;
    MMatrix driverMatrix;
    MMatrix drivenMatrix;
    MMatrix drivenInverseMatrix;
    float scale;
} TASKDATA, *PTASKDATA;
 

typedef struct ThreadData
{
    unsigned int start;
    unsigned int end;
    unsigned int numTasks;
    PTASKDATA pData;
} THREADDATA, *PTHREADDATA;


class cvWrap : public MPxDeformerNode
{
public:
						    cvWrap();
	virtual				    ~cvWrap(); 
    virtual MStatus		    deform( MDataBlock& data, MItGeometry& iter, const MMatrix& mat, unsigned int mIndex );

            MStatus         getBindInfo( MDataBlock& data );
    static void quickSort( int low, int high, std::vector<Sortable>& values );

	static  void*		    creator();
	static  MStatus		    initialize();
            void            createThreadData( int numTasks, TASKDATA *pTaskData, PTHREADDATA &pThreadData );
    static  void            createTasks( void *data, MThreadRootTask *pRoot );
    static  MThreadRetVal   threadEvaluate( void *pParam );
    static  void createMatrix( MMatrix& matrix, MPoint& origin, MVector& normal, MVector& up );
    
    static  MObject         aBindDriverGeo;
    static  MObject         aDriverGeo;
    static  MObject         aSampleVerts;
    static  MObject         aRadius;
    static  MObject         aBindInfo;
    static  MObject         aSampleComponents;
    static  MObject         aSampleWeights;
    static  MObject         aBindMatrix;
    static  MObject         aTriangleVerts;
    static  MObject         aBarycentricWeights;
    static  MObject         aNumTasks;
    static  MObject         aDirty;
    static  MObject         aScale;
    static  MObject         aCVsInU;
    static  MObject         aCVsInV;
    static  MObject         aFormInU;
    static  MObject         aFormInV;
    static  MObject         aDegreeV;



	static	MTypeId		    id;

private:
    std::vector< std::vector <Sortable> > m_sortedWeights;
    MMatrixArray                m_bindMatrices;
    std::vector<MIntArray>      m_triangleVerts;
    std::vector<MFloatArray>    m_barycentricWeights;
    TASKDATA                    m_taskData;

};

#endif
