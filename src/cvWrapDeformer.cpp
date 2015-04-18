#include "cvWrapDeformer.h"

MTypeId     cvWrap::id( 0x0011580B );

MObject     cvWrap::aBindDriverGeo;
MObject     cvWrap::aDriverGeo;
MObject     cvWrap::aSampleVerts;
MObject     cvWrap::aRadius;
MObject     cvWrap::aBindInfo;
MObject     cvWrap::aSampleComponents;
MObject     cvWrap::aSampleWeights;
MObject     cvWrap::aBindMatrix;
MObject     cvWrap::aTriangleVerts;
MObject     cvWrap::aBarycentricWeights;
MObject     cvWrap::aNumTasks;
MObject     cvWrap::aDirty;
MObject     cvWrap::aScale;
MObject     cvWrap::aCVsInU;
MObject     cvWrap::aCVsInV;
MObject     cvWrap::aFormInU;
MObject     cvWrap::aFormInV;
MObject     cvWrap::aDegreeV;


cvWrap::cvWrap() 
{
    MThreadPool::init();
}

cvWrap::~cvWrap() 
{
    MThreadPool::release();
    for ( unsigned int i = 0; i < m_sortedWeights.size(); i++ )
    {
        m_sortedWeights[i].clear();
    }
    m_sortedWeights.clear();
    m_triangleVerts.clear();
    m_barycentricWeights.clear();
}


void* cvWrap::creator() { return new cvWrap(); }


MStatus cvWrap::deform( MDataBlock& data,
        MItGeometry& itGeo,
        const MMatrix& localToWorldMatrix,
        unsigned int geomIndex )
{
    MStatus status;
    // Get envelope
	m_taskData.envelope = data.inputValue( envelope ).asFloat();
    int numTasks = data.inputValue( aNumTasks ).asInt();
    if ( m_taskData.envelope == 0.0f || numTasks <= 0 )
    {
        return MS::kSuccess;
    }

    // Get driver geo
    MDataHandle hDriverGeo = data.inputValue( aDriverGeo, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MObject oDriverGeo = hDriverGeo.asMesh();
    CHECK_MSTATUS_AND_RETURN_IT( status );
    if ( oDriverGeo.isNull() )
    {
        return MS::kSuccess;
    }

    bool dirty = data.inputValue( aDirty ).asBool();

    // Sorted weights are used to build a correct coordinate space.
    // There are cases when the up vector may be 0 length, in which
    // case we need to take out vertex influence starting at the lowest
    // weight.
    if ( m_taskData.sortedWeights.size() == 0 || dirty )
    {
        MDataHandle hDirty = data.outputValue( aDirty, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        hDirty.setBool( false );
        hDirty.setClean();

        status = getBindInfo( data );
        if ( status == MS::kNotImplemented )
        {
            return MS::kSuccess;
        }
        else if ( MFAIL( status ) )
        {
            CHECK_MSTATUS_AND_RETURN_IT( status );
        }
    }
    

    // Get driver geo information
    MFnMesh fnDriver( oDriverGeo, &status );
    MFloatVectorArray driverNormalsFloat;
    CHECK_MSTATUS_AND_RETURN_IT( status );
    status = fnDriver.getPoints( m_taskData.driverPoints, MSpace::kWorld );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    unsigned int numDriverPoints = m_taskData.driverPoints.length();

    // Put desired points and normals into world space
    bool* pAlreadyCalculated = new bool[numDriverPoints];
    for ( unsigned int i = 0; i < numDriverPoints; i++ )
    {
        pAlreadyCalculated[i] = false;
    }

    int index = 0;
    m_taskData.driverNormals.setLength( numDriverPoints );
    for( unsigned int i = 0; i < m_taskData.sortedWeights.size(); i++ )
    {
        for( unsigned int j = 0; j < m_taskData.sortedWeights[i].size(); j++ )
        {
            index = m_taskData.sortedWeights[i][j].index;
            if ( !pAlreadyCalculated[index] )
            {
                status = fnDriver.getVertexNormal( index, false, m_taskData.driverNormals[index], MSpace::kWorld );
                CHECK_MSTATUS_AND_RETURN_IT( status );
            }
        }
    }
    delete [] pAlreadyCalculated;
    
    m_taskData.scale = data.inputValue( aScale ).asFloat();

    m_taskData.membership.setLength( itGeo.count() );
    status = itGeo.allPositions( m_taskData.points );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    m_taskData.weights.setLength( itGeo.count() );
    int i = 0;
    
    bool isNurbsSurface = false;
    for ( itGeo.reset(); !itGeo.isDone(); itGeo.next(), i++ )
    {
		m_taskData.membership[i] = itGeo.index();
        if ( m_taskData.membership[i] >= m_sortedWeights.size() )
        {
            isNurbsSurface = true;
        }
        m_taskData.weights[i] = weightValue( data, geomIndex, itGeo.index() );
    }

    // Update membership indices for nurbs surfaces 
    // because indices from itGeo return more than cv indices.
    if ( isNurbsSurface )
    {
        int uLimit = data.inputValue( aCVsInU ).asInt();
        int cvsInV = data.inputValue( aCVsInV ).asInt();
        int vLimit = cvsInV;
        if ( !(data.inputValue( aFormInU ).asInt() == 1 &&data.inputValue( aFormInV ).asInt() == 1) )
        {
            vLimit -= data.inputValue( aDegreeV ).asInt();
        }
        int memberIndex = 0;
        int cvIndex = 0;
        for ( int u = 0; u < uLimit; u++ )
        {
            for ( int v = 0; v < vLimit; v++ )
            {
                int index = cvsInV * u + v;
                if ( index == m_taskData.membership[memberIndex] )
                {
                    m_taskData.membership[memberIndex] = cvIndex;
                    memberIndex++;
                }
                cvIndex++;
            }
        }
    }

    if ( itGeo.count() > 80000 )
    {
        MGlobal::displayError( "Demo version has 80000 vertex limit." );
        return MS::kSuccess;
    }

    

    m_taskData.drivenMatrix = localToWorldMatrix;
    m_taskData.drivenInverseMatrix = localToWorldMatrix.inverse();
    
    PTHREADDATA pThreadData = NULL;
    createThreadData( numTasks, &m_taskData, pThreadData );
    MThreadPool::newParallelRegion( createTasks, (void *)pThreadData );

    status = itGeo.setAllPositions( m_taskData.points );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    if ( pThreadData )
    {
        delete [] pThreadData;
    }
    return MS::kSuccess;
}


MStatus cvWrap::getBindInfo( MDataBlock& data )
{
    MStatus status;

    MArrayDataHandle hSampleWeights = data.inputArrayValue( aSampleWeights );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    unsigned int numVerts = hSampleWeights.elementCount();
    if ( numVerts == 0 )
    {
        return MS::kNotImplemented;
    }

    MArrayDataHandle hComponents = data.inputArrayValue( aSampleComponents, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MArrayDataHandle hBindMatrix = data.inputArrayValue( aBindMatrix, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MArrayDataHandle hTriangleVerts = data.inputArrayValue( aTriangleVerts, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MArrayDataHandle hBarycentricWeights = data.inputArrayValue( aBarycentricWeights, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MObject oWeightData, oComponentData, oNumericData;
    MIntArray compIds, tempIds;
    MFnSingleIndexedComponent fnSingleComp;
    MFnComponentListData fnCompData;
    MFnNumericData fnNumericData;
    m_taskData.sortedWeights.resize( numVerts );
    m_taskData.bindMatrices.setLength( numVerts );
    m_taskData.triangleVerts.resize( numVerts );
    m_taskData.barycentricWeights.resize( numVerts );
    MDoubleArray weights, normalizedWeights;
    double totalWeight = 0.0;
    for ( unsigned int i = 0; i < numVerts; i++ )
    {
        // Get component ids
        status = hComponents.jumpToArrayElement( i );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        oComponentData = hComponents.inputValue().data();
        status = fnCompData.setObject( oComponentData );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        compIds.clear();
        for ( unsigned int j = 0; j < fnCompData.length(); j++ )
        {
            if ( fnCompData[j].hasFn( MFn::kSingleIndexedComponent ) )
            {
                status = fnSingleComp.setObject( fnCompData[j] );
                CHECK_MSTATUS_AND_RETURN_IT( status );
                fnSingleComp.getElements( tempIds );
                for( unsigned int k = 0; k < tempIds.length(); k++ )
                {
                    compIds.append( tempIds[k] );
                }
            }
        }
        m_taskData.sortedWeights[i].resize( compIds.length() );
        for ( unsigned int j = 0; j < compIds.length(); j++ )
        {
            m_taskData.sortedWeights[i][j].index = compIds[j];
        }

        // Get sample weights
        status = hSampleWeights.jumpToArrayElement( i );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        oWeightData = hSampleWeights.inputValue().data();
        MFnDoubleArrayData fnDoubleData( oWeightData, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        weights = fnDoubleData.array();
        totalWeight = 0.0;
        for ( unsigned int j = 0; j < weights.length(); j++ )
        {
            totalWeight += weights[j];
            m_taskData.sortedWeights[i][j].weight = weights[j];

        }

        // Calculate normalized weights
        normalizedWeights.setLength( weights.length() );
        if ( totalWeight == 0.0 )
        {
            CHECK_MSTATUS_AND_RETURN_IT( MS::kFailure );
        }
        for ( unsigned int j = 0; j < weights.length(); j++ )
        {
            normalizedWeights[j] = weights[j] / totalWeight;
            m_taskData.sortedWeights[i][j].normalizedWeight = normalizedWeights[j];
        }

        // Sort weights lowest to highest so we can take the low weight influences away first if
        // the up vector needs adjustment
        quickSort( 0, m_taskData.sortedWeights[i].size() - 1, m_taskData.sortedWeights[i] );

        // Get bind matrix
        status = hBindMatrix.jumpToArrayElement( i );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        m_taskData.bindMatrices[i] = hBindMatrix.inputValue().asMatrix();

        // Get triangle vertex binding
        status = hTriangleVerts.jumpToArrayElement( i );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        int3& verts = hTriangleVerts.inputValue( &status ).asInt3();
        CHECK_MSTATUS_AND_RETURN_IT( status );
        m_taskData.triangleVerts[i].setLength( 3 );
        m_taskData.triangleVerts[i][0] = verts[0];
        m_taskData.triangleVerts[i][1] = verts[1];
        m_taskData.triangleVerts[i][2] = verts[2];

        // Get barycentric weights
        status = hBarycentricWeights.jumpToArrayElement( i );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        float3& baryWeights = hBarycentricWeights.inputValue( &status ).asFloat3();
        CHECK_MSTATUS_AND_RETURN_IT( status );
        m_taskData.barycentricWeights[i].setLength( 3 );
        m_taskData.barycentricWeights[i][0] = baryWeights[0];
        m_taskData.barycentricWeights[i][1] = baryWeights[1];
        m_taskData.barycentricWeights[i][2] = baryWeights[2];
    }
    return MS::kSuccess;
}


void cvWrap::quickSort( int low, int high, std::vector<Sortable>& values )
{
    int i = low;
    int j = high;
    double pivot = values[(low + high) / 2].weight;
    Sortable temp;
    while ( i <= j )
    {
        while ( values[i].weight < pivot )
        {
            if ( i + 1 >= values.size() )
            {
                break;
            }
            i++;
        }
        while ( values[j].weight > pivot )
        {
            if ( j <= 0 )
            {
                break;
            }
            j--;
        }
        if ( i <= j )
        {
            temp = values[i];
            values[i] = values[j];
            values[j] = temp;
            i++;
            j--;
        }
    }
    if ( low < j )
    {
        quickSort( low, j, values );
    }
    if ( i < high )
    {
        quickSort( i, high, values );
    }
}


void cvWrap::createThreadData( int numTasks, TASKDATA *pTaskData, PTHREADDATA &pThreadData )
{
    if ( pThreadData )
    {
        delete [] pThreadData;
    }
    pThreadData = new THREADDATA[numTasks];
    unsigned int taskLength = (pTaskData->points.length() + numTasks - 1) / numTasks;
    unsigned int start = 0;
    unsigned int end = taskLength;

    int lastTask = numTasks - 1;
	for( int i = 0; i < numTasks; i++ )
	{
        if ( i == lastTask )
        {
            end = pTaskData->points.length();
        }
        pThreadData[i].start = start;
        pThreadData[i].end = end;
        pThreadData[i].numTasks = numTasks;
        pThreadData[i].pData = pTaskData;

        start += taskLength;
        end += taskLength;
	}
}


void cvWrap::createTasks( void *pData, MThreadRootTask *pRoot )
{
    THREADDATA * pThreadData = (THREADDATA *)pData;

	if ( pThreadData )
	{
        int numTasks = pThreadData->numTasks;
		for( int i = 0; i < numTasks; i++ )
		{
            MThreadPool::createTask( threadEvaluate, (void *)&pThreadData[i], pRoot );
		}
		MThreadPool::executeAndJoin( pRoot );
	}
}


MThreadRetVal cvWrap::threadEvaluate( void *pParam )
{
    PTHREADDATA pThreadData = (PTHREADDATA)(pParam);
    PTASKDATA pData = pThreadData->pData;
    float env = pThreadData->pData->envelope;
    MIntArray& membership = pData->membership;
    MFloatArray& weights = pData->weights;
    MPointArray& points = pData->points;
    unsigned int numDeformVerts = points.length();

    MPointArray& driverPoints = pData->driverPoints;
    MVectorArray& driverNormals = pData->driverNormals;
    std::vector < std::vector <Sortable> >& sortedWeights = pData->sortedWeights;
    MMatrixArray& bindMatrices = pData->bindMatrices;
    std::vector <MIntArray>& triangleVerts = pData->triangleVerts;
    std::vector <MFloatArray>& barycentricWeights = pData->barycentricWeights;
    MMatrix& driverMatrix = pData->driverMatrix;
    MMatrix& drivenMatrix = pData->drivenMatrix;
    MMatrix& drivenInverseMatrix = pData->drivenInverseMatrix;
    float scale = pThreadData->pData->scale;

    unsigned int taskStart = pThreadData->start;
    unsigned int taskEnd = pThreadData->end;
    int index;

    unsigned int numSamples;
    MPoint origin, newPt, pt;
    MVector normal, up;
    MMatrix matrix;
	for ( unsigned int i = taskStart; i < taskEnd; i++ )
    {
        if ( i >= numDeformVerts )
        {
            continue;
        }
        index = membership[i];
        numSamples = sortedWeights[index].size();

        // Reconstruct origin
        origin = (MVector( driverPoints[triangleVerts[index][0]] ) * barycentricWeights[index][0]) + 
            (MVector( driverPoints[triangleVerts[index][1]] ) * barycentricWeights[index][1]) +
            (MVector( driverPoints[triangleVerts[index][2]] ) * barycentricWeights[index][2]);


        // Reconstruct normal and up vector
        normal = MVector::zero;
        up = MVector::zero;
        for ( unsigned int j = 0; j < numSamples; j++ )
        {
            normal += driverNormals[sortedWeights[index][j].index] * sortedWeights[index][j].normalizedWeight;
            up += (driverPoints[sortedWeights[index][j].index] - origin) * sortedWeights[index][j].normalizedWeight;
        }
        normal.normalize();
        normal *= scale;

        // If the up and normal or parallel or the up has 0 length, take out some influence.
        if ( up * normal == 1.0 || up.length() < 0.0001 )
        {
            for ( unsigned int j = 0; j < numSamples; j++ )
            {
                if ( up * normal != 1.0 && up.length() > 0.0001 )
                {
                    break;
                }
                up -= (driverPoints[sortedWeights[index][j].index] - origin) * sortedWeights[index][j].normalizedWeight;
            }
        }   

        up.normalize();
        
        createMatrix( matrix, origin, normal, up );
        
        pt = points[i];
        newPt = ((pt  * drivenMatrix) * (bindMatrices[index] * matrix)) * drivenInverseMatrix;
        points[i] = pt + ((newPt - pt) * weights[i] * env);
    }
    return 0;
}


void cvWrap::createMatrix( MMatrix& matrix, MPoint& origin, MVector& normal, MVector& up )
{
    MVector x = normal ^ up;
    up = normal ^ x;
    matrix[0][0] = x.x;         matrix[0][1] = x.y;         matrix[0][2] = x.z;      matrix[0][3] = 0.0;
    matrix[1][0] = normal.x;    matrix[1][1] = normal.y;    matrix[1][2] = normal.z; matrix[1][3] = 0.0;
    matrix[2][0] = up.x;        matrix[2][1] = up.y;        matrix[2][2] = up.z;     matrix[2][3] = 0.0;
    matrix[3][0] = origin.x;    matrix[3][1] = origin.y;    matrix[3][2] = origin.z; matrix[3][3] = 1.0;
}


MStatus cvWrap::initialize()
{
    MFnCompoundAttribute    cAttr;
    MFnMatrixAttribute      mAttr;
    MFnMessageAttribute     meAttr;
    MFnTypedAttribute       tAttr;
    MFnGenericAttribute     gAttr;
    MFnNumericAttribute     nAttr;
	MStatus				    status;

    aDriverGeo = tAttr.create( "driver", "driver", MFnData::kMesh );
    addAttribute( aDriverGeo );
    attributeAffects( aDriverGeo, outputGeom );

    aBindDriverGeo = meAttr.create( "bindMesh", "bindMesh" );
    addAttribute( aBindDriverGeo );

    aSampleComponents = tAttr.create( "sampleComponents", "sampleComponents", MFnData::kComponentList );
    tAttr.setArray( true );
    addAttribute( aSampleComponents );

    aSampleWeights = tAttr.create( "sampleWeights", "sampleWeights", MFnData::kDoubleArray );
    tAttr.setArray( true );
    addAttribute( aSampleWeights );

    aBindMatrix = mAttr.create( "bindMatrix", "bindMatrix" );
    mAttr.setDefault( MMatrix::identity );
    mAttr.setArray( true );
    addAttribute( aBindMatrix );

    aTriangleVerts = nAttr.create( "triangleVerts", "triangleVerts", MFnNumericData::k3Int );
    nAttr.setArray( true );
    addAttribute( aTriangleVerts );

    aBarycentricWeights = nAttr.create( "barycentricWeights", "barycentricWeights", MFnNumericData::k3Float );
    nAttr.setArray( true );
    addAttribute( aBarycentricWeights );

    aScale = nAttr.create( "scale", "scale", MFnNumericData::kFloat, 1.0 );
    nAttr.setKeyable( true );
    addAttribute( aScale );
    attributeAffects( aScale, outputGeom );

    aNumTasks = nAttr.create( "numTasks", "numTasks", MFnNumericData::kInt, 32 );
    nAttr.setMin( 1 );
    addAttribute( aNumTasks );

    aDirty = nAttr.create( "dirty", "dirty", MFnNumericData::kBoolean );
    addAttribute( aDirty );
    attributeAffects( aDirty, outputGeom );

    aCVsInU = nAttr.create( "cvsInU", "cvsInU", MFnNumericData::kInt, 0 );
    nAttr.setHidden( true );
    addAttribute( aCVsInU );

    aCVsInV = nAttr.create( "cvsInV", "cvsInV", MFnNumericData::kInt, 0 );
    nAttr.setHidden( true );
    addAttribute( aCVsInV );

    aFormInU = nAttr.create( "formInU", "formInU", MFnNumericData::kInt, 0 );
    nAttr.setHidden( true );
    addAttribute( aFormInU );

    aFormInV = nAttr.create( "formInV", "formInV", MFnNumericData::kInt, 0 );
    nAttr.setHidden( true );
    addAttribute( aFormInV );

    aDegreeV = nAttr.create( "degreesV", "degreesV", MFnNumericData::kInt, 0 );
    nAttr.setHidden( true );
    addAttribute( aDegreeV );

	MGlobal::executeCommand( "makePaintable -attrType multiFloat -sm deformer cvWrap weights" );
    
    return MS::kSuccess;
}


