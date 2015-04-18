
#include "cvWrapCmd.h"
#include "cvWrapRebindCmd.h"
#include "cvWrapDeformer.h"


cvWrapRebindCmd::cvWrapRebindCmd()
{
}


cvWrapRebindCmd::~cvWrapRebindCmd()
{
    m_origSampleVerts.clear();
    m_origSampleWeights.clear();
    m_origTriangleVertices.clear();
    m_origBarycentricCoordinates.clear();
}


MSyntax cvWrapRebindCmd::newSyntax()
{
    MSyntax syntax;
    syntax.addFlag( "-r", "-radius", MSyntax::kDouble );
    syntax.addFlag( "-wn", "-wrapNode", MSyntax::kString );
    syntax.addFlag( "-wmm", "-weightMapMultiply" );

    syntax.setObjectType( MSyntax::kSelectionList, 1, 255 );
    syntax.useSelectionAsDefault( true );
    return syntax;
}


void* cvWrapRebindCmd::creator()
{
    return new cvWrapRebindCmd;
}


bool cvWrapRebindCmd::isUndoable() const
{
   return true;
}


MStatus cvWrapRebindCmd::doIt( const MArgList& args )
{
    MStatus status;
    m_radius = 0.1;
    MArgDatabase argData( syntax(), args );
    if ( argData.isFlagSet( "-r" ) )
    {
        m_radius = argData.flagArgumentDouble( "-r", 0, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
    }
    if ( argData.isFlagSet( "-wn" ) )
    {
        m_wrapNode = argData.flagArgumentDouble( "-wn", 0, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
    }
    m_weightMapMultiply = argData.isFlagSet( "-wmm" );
    argData.getObjects( m_selectionList );

    // Get selection paths
    status = getPathsToDriverAndDriven( m_pathDriver, m_oCompDriver,
        m_pathDriven, m_oCompDriven );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    // Get vertices to rebind
    MFnSingleIndexedComponent fnDrivenComp( m_oCompDriven, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    status = fnDrivenComp.getElements( m_selectedVerts );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    // Get selected driver faces
    MFnSingleIndexedComponent fnDriverComp( m_oCompDriver, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    status = fnDriverComp.getElements( m_selectedFaces );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    m_origSampleVerts.resize( m_selectedVerts.length() );
    m_origSampleWeights.resize( m_selectedVerts.length() );
    m_origTriangleVertices.resize( m_selectedVerts.length() );
    m_origBarycentricCoordinates.resize( m_selectedVerts.length() );
    m_origBindMatrices.setLength( m_selectedVerts.length() );

    status = redoIt();
    CHECK_MSTATUS_AND_RETURN_IT( status );

    return MS::kSuccess;
}


MStatus cvWrapRebindCmd::redoIt()
{
    MStatus status;

    // Get wrap node
    MObject oWrapNode;

    if ( m_wrapNode.length() == 0 )
    {
        status = getWrapNode( oWrapNode );
        if ( MFAIL( status ) )
        {
            MGlobal::displayError( "Unable to find wrap node.  Make sure both wrapped and wrapper meshes are selected." );
        }
        CHECK_MSTATUS_AND_RETURN_IT( status );
        MFnDependencyNode fnWrapNode( oWrapNode, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        m_wrapNode = fnWrapNode.name();
    }
    else
    {
        status = getMObject( m_wrapNode, oWrapNode );
        CHECK_MSTATUS_AND_RETURN_IT( status );
    }


    if ( m_weightMapMultiply )
    {
        //status = multiplyWeightsWithMap( pathDriven, oCompDriven, oWrapNode );
    }
    else
    {
        status = rebindVertices( oWrapNode );
        CHECK_MSTATUS_AND_RETURN_IT( status );
    }

    MPlug plugDirty( oWrapNode, cvWrap::aDirty );
    plugDirty.setBool( true );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    MGlobal::setActiveSelectionList( m_selectionList );
    

    return MS::kSuccess;
}


MStatus cvWrapRebindCmd::getPathsToDriverAndDriven(
        MDagPath& pathDriver, MObject& oCompDriver,
        MDagPath& pathDriven, MObject& oCompDriven )
{
    MStatus status;
    // Get selection paths
    MItSelectionList itSel( m_selectionList, MFn::kInvalid, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    for ( ; !itSel.isDone(); itSel.next() )
    {
        if ( !pathDriven.isValid() )
        {
            status = itSel.getDagPath( pathDriven, oCompDriven );
            CHECK_MSTATUS_AND_RETURN_IT( status );
        }
        else
        {
            status = itSel.getDagPath( pathDriver, oCompDriver );
            CHECK_MSTATUS_AND_RETURN_IT( status );
        }
    }

    // If the transform was selected, get the shape.
    if ( pathDriver.node().hasFn( MFn::kTransform ) )
    {
        unsigned int childCount = pathDriver.childCount();
        for ( unsigned int j = 0; j < childCount; j++ )
        {
            MObject oChild = pathDriver.child( j, &status );
            CHECK_MSTATUS_AND_RETURN_IT( status );
            pathDriver.push( oChild );

            MFnDagNode fnDag( pathDriver );
            if ( fnDag.isIntermediateObject() )
            {
                pathDriver.pop();
                continue;
            }

            if ( pathDriver.node().apiType() == MFn::kMesh 
                    || pathDriver.node().apiType() == MFn::kNurbsSurface 
                    || pathDriver.node().apiType() == MFn::kNurbsCurve )
            {
                break;
            }
        }
    }
    return MS::kSuccess;
}


MStatus cvWrapRebindCmd::getWrapNode( MObject& oWrapNode )
{
    MStatus status;
    MFnDagNode fnDriven( m_pathDriven, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MPlug plugShape = fnDriven.findPlug( "inMesh", false, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    MItDependencyGraph itDep( plugShape, MFn::kInvalid, MItDependencyGraph::kUpstream, 
        MItDependencyGraph::kDepthFirst, MItDependencyGraph::kPlugLevel, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MFnDependencyNode fnNode;
    for ( ; !itDep.isDone(); itDep.next() )
    {
        // The node we want has wrapMesh.outMesh connected to node.driver
        fnNode.setObject( itDep.currentItem() );
        if ( fnNode.typeId() != cvWrap::id )
        {
            // Not a cvWrap
            continue;
        }

        MPlug plugDriver = fnNode.findPlug( "driver", false, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Check if wrapMesh is connected to this plug
        MPlugArray plugArray;
        plugDriver.connectedTo( plugArray, true, false, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        if ( plugArray.length() == 0 )
        {
            // No plugs connected
            continue;
        }
        if ( !plugArray[0].node().hasFn( MFn::kDagNode ) )
        {
            // Wrap mesh may be a direct connection from another deformer
            continue;
        }
        MFnDagNode fnDag( plugArray[0].node(), &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        MDagPath pathMesh;
        fnDag.getPath( pathMesh );
        if ( !( pathMesh == m_pathDriver ) )
        {
            // No match
            continue;
        }
        // If we get here, we found a match
        oWrapNode = itDep.currentItem();
        return MS::kSuccess;
    }
    return MS::kFailure;
}


//MStatus cvWrapRebindCmd::multiplyWeightsWithMap( MDagPath& pathDriven, MObject& oComponents, MObject& oWrapNode )
//{
//    MStatus status;
//    // Get current color set of driver
//    MPlug plugInDriverMesh( oWrapNode, cvWrap::aDriverGeo );
//    MPlug plugDriverMesh;
//    status = imdUtility::GetConnectedPlug( plugInDriverMesh, plugDriverMesh );
//    CHECK_MSTATUS_AND_RETURN_IT( status );
//    MObject oDriverMesh = plugDriverMesh.node();
//    MFnMesh fnDriverMesh( oDriverMesh, &status );
//    CHECK_MSTATUS_AND_RETURN_IT( status );
//    MString currentColorSetName = fnDriverMesh.currentColorSetName( 0 );
//    if ( currentColorSetName.numChars() == 0 )
//    {
//        return MS::kSuccess;
//    }
//    MColor defaultColor( 1.0f, 1.0f, 1.0f, 1.0f );
//    MColorArray colors;
//    status = fnDriverMesh.getVertexColors( colors, &currentColorSetName, &defaultColor );
//    CHECK_MSTATUS_AND_RETURN_IT( status );
//
//    MPlug plugSampleWeights( oWrapNode, cvWrap::aSampleWeights );
//    MPlug plugSampleVerts( oWrapNode, cvWrap::aSampleComponents );
//
//
//    m_origSampleWeights.resize( m_selectedVerts.length() );
//    MPlug plugTemp;
//    MObject oTemp;
//    MDoubleArray tempWeights;
//    MFnSingleIndexedComponent fnSampleVertComp;
//    MFnComponentListData fnCompData;
//    MIntArray sampleVertIds;
//    for ( unsigned int i = 0; i < m_selectedVerts.length(); i++ )
//    {
//        // Get sample vert ids.
//        plugTemp = plugSampleVerts.elementByPhysicalIndex( m_selectedVerts[i], &status );
//        CHECK_MSTATUS_AND_RETURN_IT( status );
//        oTemp = plugTemp.asMObject();
//        status = fnCompData.setObject( oTemp );
//        CHECK_MSTATUS_AND_RETURN_IT( status );
//        status = fnSampleVertComp.setObject( fnCompData[0] );
//        CHECK_MSTATUS_AND_RETURN_IT( status );
//        sampleVertIds.clear();
//        status = fnSampleVertComp.getElements( sampleVertIds );
//        CHECK_MSTATUS_AND_RETURN_IT( status );
//
//        // Get current sample weights
//        plugTemp = plugSampleWeights.elementByPhysicalIndex( m_selectedVerts[i], &status );
//        CHECK_MSTATUS_AND_RETURN_IT( status );
//        oTemp = plugTemp.asMObject();
//        MFnDoubleArrayData fnDoubleData( oTemp );
//        tempWeights = fnDoubleData.array();
//        double * w = new double[tempWeights.length()];
//        tempWeights.get( w );
//        m_origSampleWeights[i] = MDoubleArray( w, tempWeights.length() );
//        SAFE_DELETE_ARRAY( w );
//
//        // Multiply each sample weight
//        if ( sampleVertIds.length() == tempWeights.length() )
//        {
//            for ( unsigned int j = 0; j < tempWeights.length(); j++ )
//            {
//                tempWeights[j] *= colors[sampleVertIds[j]].r;
//            }
//        }
//        else
//        {
//            cerr << "Sample vertex Ids != sample weight values.\n";
//        }
//
//        oTemp = fnDoubleData.create( tempWeights, &status );
//        CHECK_MSTATUS_AND_RETURN_IT( status );
//        plugTemp.setMObject( oTemp );
//    }
//    return MS::kSuccess;
//}


MStatus cvWrapRebindCmd::rebindVertices( MObject& oWrapNode )
{
    MStatus status;

    // Get the bind mesh
    MDagPath pathBindMesh;
    status = getBindMesh( oWrapNode, pathBindMesh );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MFnMesh fnBindMesh( pathBindMesh, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    // Get intermediate driven mesh
    MDagPath pathIntermediate;
    status = getIntermediateObject( m_pathDriven, pathIntermediate );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    // Create driver mesh based off of specified faces
    MDagPath pathDriverSubset;
    status = createDriverSubsetMesh( pathBindMesh, pathDriverSubset );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    
    // Create meshIntersectors
    MMeshIntersector subsetIntersector;
    MObject oDriverSubset = pathDriverSubset.node();
    MMatrix driverSubsetMatrix = pathDriverSubset.inclusiveMatrix();
    status = subsetIntersector.create( oDriverSubset, driverSubsetMatrix );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    MMeshIntersector intersector;
    MMatrix driverBindMatrix = pathBindMesh.inclusiveMatrix();
    MObject oBindMesh = pathBindMesh.node();
    status = intersector.create( oBindMesh, driverBindMatrix );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    MFnMesh fnDriven( m_pathDriven, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MFnMesh fnDrivenOrig( pathIntermediate, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MPointArray origPoints, bindPoints;
    status = fnDrivenOrig.getPoints( origPoints, MSpace::kWorld );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    status = fnBindMesh.getPoints( bindPoints, MSpace::kWorld );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MFloatVectorArray normals;
    status = fnBindMesh.getVertexNormals( false, normals, MSpace::kWorld );
    CHECK_MSTATUS_AND_RETURN_IT( status );


    std::vector<std::set<int> > adjacency;
    status = CVWrapCmd::getAdjacency( m_pathDriver, adjacency );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    MPoint pt, closestPoint;
    MPointOnMesh pointOnMesh;
    int triangleVertices[3];
    float a, b, c;
    double totalWeight, radius;
    MIntArray vertexList, sampleIds;
    int faceId, triangleId;
    MFnSingleIndexedComponent fnComp;
    MFnComponentListData fnCompData;
    MSelectionList selectionList;
    MDagPath pathTemp;
    MObject oComp, oCompData, oNumericData, oMatrixData, oDoubleData;
    MDoubleArray weights, normalizedWeights;
    MPoint origin, localPoint;
    MVector normal, up;
    MPlug plugSampleWeights( oWrapNode, cvWrap::aSampleWeights );
    MPlug plugSampleVerts( oWrapNode, cvWrap::aSampleComponents );
    MPlug plugSampleBindMatrix( oWrapNode, cvWrap::aBindMatrix );
    MPlug plugTriangleVerts( oWrapNode, cvWrap::aTriangleVerts );
    MPlug plugBarycentricWeights( oWrapNode, cvWrap::aBarycentricWeights );
    MPlug plugTemp;

    MFnSingleIndexedComponent fnSampleVertComp;
    MFnNumericData fnNumericData;
    MFnMatrixData fnMatrixData;
    MMatrix matrix;
    MObject oTemp;
    MDoubleArray tempWeights;

    bool progressbar = MGlobal::mayaState() == MGlobal::kInteractive;
    // Start progress bar
    StartProgress("Rebinding wrap...", m_selectedVerts.length());
    for (unsigned int i = 0; i < m_selectedVerts.length(); i++) {
      if (ProgressCancelled()) {
        break;
      }

        // Closest point to subset
        pt = origPoints[m_selectedVerts[i]];
        status = subsetIntersector.getClosestPoint( pt, pointOnMesh );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        closestPoint = MPoint( pointOnMesh.getPoint() ) * driverSubsetMatrix;

        // Closest point to bind mesh from subset
        status = intersector.getClosestPoint( closestPoint, pointOnMesh );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        closestPoint = MPoint( pointOnMesh.getPoint() ) * driverBindMatrix;
        faceId = pointOnMesh.faceIndex();
        triangleId = pointOnMesh.triangleIndex();

        // Get vertices of closest triangle
        vertexList.clear();
        status = fnBindMesh.getPolygonTriangleVertices( faceId, triangleId, triangleVertices );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        //oComp = fnComp.create( MFn::kMeshVertComponent, &status );
        //CHECK_MSTATUS_AND_RETURN_IT( status );
        //status = fnComp.addElements( vertexList );
        //CHECK_MSTATUS_AND_RETURN_IT( status );

        // Get barycentric coordinates of closestPoint
        CVWrapCmd::getBarycentricCoordinates( closestPoint,
                bindPoints[triangleVertices[0]],
                bindPoints[triangleVertices[1]],
                bindPoints[triangleVertices[2]], a, b, c );

        // Get vertices of closest face
        vertexList.clear();
        status = fnBindMesh.getPolygonVertices( faceId, vertexList );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Crawl surface to get sample weights
        std::map<int, double> distances;
        radius = m_radius;
        do
        {
            distances.clear();
            for(unsigned int j = 0; j < vertexList.length(); ++j) {
                status = CVWrapCmd::crawlSurface( vertexList[j], bindPoints, distances,
                        0.0, closestPoint, radius, adjacency );
                CHECK_MSTATUS_AND_RETURN_IT( status );
            }

            // Get sample weights
            CVWrapCmd::calculateSampleWeights( distances, sampleIds, weights, normalizedWeights,
                    closestPoint, bindPoints, totalWeight );
            radius *= 2.0;
        } while ( distances.size() == 0 || totalWeight == 0.0 );

        // Calculate origin
        origin = closestPoint;

        // Calculate normal and up
        normal = MVector::zero;
        up = MVector::zero;
        for ( unsigned int j = 0; j < weights.length(); j++ )
        {
            normal += MVector( normals[sampleIds[j]] ) * normalizedWeights[j];
            up += (bindPoints[sampleIds[j]] - origin) * normalizedWeights[j];
        }

        // Adjust up if it's parallel to normal or if it's zero length
        if ( up * normal == 1.0 || up.length() < 0.0001 )
        {
            // Sort weights low to high.
            std::vector<Sortable> sortedWeights;
            sortedWeights.resize( weights.length() );
            for ( unsigned int i = 0; i < weights.length(); i++ )
            {
                sortedWeights[i].weight = weights[i];
                sortedWeights[i].index = sampleIds[i];
                sortedWeights[i].normalizedWeight = normalizedWeights[i];
            }
            cvWrap::quickSort( 0, weights.length() - 1, sortedWeights );

            for ( unsigned int i = 0; i < weights.length(); i++ )
            {
                if ( up * normal != 1.0 && up.length() > 0.0001 )
                {
                    break;
                }
                up -= (bindPoints[sortedWeights[i].index] - origin) * sortedWeights[i].normalizedWeight;
            }
        }
        up.normalize();
        normal.normalize();
        
        // Store sample vert ids.
        plugTemp = plugSampleVerts.elementByPhysicalIndex( m_selectedVerts[i], &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Get the current samples ids so we can restore them on undo
        oTemp = plugTemp.asMObject();
        status = fnCompData.setObject( oTemp );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = fnSampleVertComp.setObject( fnCompData[0] );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = fnSampleVertComp.getElements( m_origSampleVerts[i] );

        // Store the new sample ids
        oComp = fnComp.create( MFn::kMeshVertComponent, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = fnComp.addElements( sampleIds );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        oCompData = fnCompData.create();
        fnCompData.add( oComp );
        plugTemp.setMObject( oCompData );

        // Store sample weights
        plugTemp = plugSampleWeights.elementByPhysicalIndex( m_selectedVerts[i], &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Get the current weights
        oTemp = plugTemp.asMObject();
        MFnDoubleArrayData fnDoubleData( oTemp );
        tempWeights = fnDoubleData.array();
        double* w = new double[tempWeights.length()];
        tempWeights.get( w );
        m_origSampleWeights[i] = MDoubleArray( w, tempWeights.length() );
        delete [] w;

        // Set the new weights
        oDoubleData = fnDoubleData.create( weights, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        plugTemp.setMObject( oDoubleData );

        // Store bind matrix
        plugTemp = plugSampleBindMatrix.elementByPhysicalIndex( m_selectedVerts[i], &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Get the current bind matrix
        oTemp = plugTemp.asMObject();
        status = fnMatrixData.setObject( oTemp );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        m_origBindMatrices[i] = fnMatrixData.matrix();

        // Store the new bind matrix
        cvWrap::createMatrix( matrix, origin, normal, up );
        matrix = matrix.inverse();
        oMatrixData = fnMatrixData.create( matrix, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        plugTemp.setMObject( oMatrixData );

        // Store triangle vertices
        plugTemp = plugTriangleVerts.elementByPhysicalIndex( m_selectedVerts[i], &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Get the current triangle vertices
        oTemp = plugTemp.asMObject();
        status = fnNumericData.setObject( oTemp );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        m_origTriangleVertices[i].setLength( 3 );
        fnNumericData.getData3Int( m_origTriangleVertices[i][0], m_origTriangleVertices[i][1], m_origTriangleVertices[i][2] );

        // Store the new triangle vertices
        oNumericData = fnNumericData.create( MFnNumericData::k3Int, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = fnNumericData.setData3Int( triangleVertices[0], triangleVertices[1], triangleVertices[2] );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        plugTemp.setMObject( oNumericData );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Store barycentric coordinates
        plugTemp = plugBarycentricWeights.elementByPhysicalIndex( m_selectedVerts[i], &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Get the current bary coords
        oTemp = plugTemp.asMObject();
        status = fnNumericData.setObject( oTemp );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        m_origBarycentricCoordinates[i].setLength( 3 );
        fnNumericData.getData3Float( m_origBarycentricCoordinates[i][0], m_origBarycentricCoordinates[i][1], m_origBarycentricCoordinates[i][2] );

        // Store the new bary coords
        oNumericData = fnNumericData.create( MFnNumericData::k3Float, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = fnNumericData.setData3Float( a, b, c );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        plugTemp.setMObject( oNumericData );
        CHECK_MSTATUS_AND_RETURN_IT( status );


        // Increment the progress bar
        StepProgress(1);
    }
    EndProgress();

    // Delete the subset mesh
    pathDriverSubset.pop();
    MFnDagNode fnTemp( pathDriverSubset );
    status = MGlobal::executeCommand( "delete " + fnTemp.partialPathName() );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    return MS::kSuccess;
}


MStatus cvWrapRebindCmd::getBindMesh( MObject& oWrapNode, MDagPath& pathBindMesh )
{
    MStatus status;

    // Get the bind mesh connected to the message attribute of the wrap deformer
    MPlug plugBindMesh( oWrapNode, cvWrap::aBindDriverGeo );
    MPlug connectedPlug;
    status = getConnectedPlug( plugBindMesh, connectedPlug );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MObject oBindMesh = connectedPlug.node();
    MDagPath::getAPathTo( oBindMesh, pathBindMesh );

    return MS::kSuccess;
}


MStatus cvWrapRebindCmd::getConnectedPlug( MPlug& plug, MPlug& connectedPlug )
{
    MStatus status;
    MPlugArray plugs;
    if ( !plug.connectedTo( plugs, true, false, &status ) )
    {
        return MS::kFailure;
    }
    CHECK_MSTATUS_AND_RETURN_IT( status );
    if ( plugs.length() == 0 )
    {
        return MS::kFailure;
    }
    connectedPlug = plugs[0];
    return MS::kSuccess;
}


MStatus cvWrapRebindCmd::getIntermediateObject( MDagPath& pathShape, MDagPath& pathIntermediate )
{
    MStatus status;
    if ( pathShape.node().hasFn( MFn::kMesh ) )
    {
        status = pathShape.pop();
    }
    if ( !pathShape.node().hasFn( MFn::kTransform ) )
    {
        CHECK_MSTATUS_AND_RETURN_IT( MS::kFailure );
    }

    unsigned int childCount = pathShape.childCount();
    MFnDagNode fnTransform( pathShape );
    MString name = fnTransform.name();
    bool hasShape = false;
    for ( unsigned int i = 0; i < childCount; i++ )
    {
        MObject oChild = fnTransform.child( i, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = pathShape.push( oChild );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        MFnDagNode fnDag( pathShape );
        if ( fnDag.isIntermediateObject() )
        {
            MPlugArray connectedPlugs;
            status = fnDag.getConnections( connectedPlugs );
            if ( status == MS::kSuccess )
            {
                if ( connectedPlugs.length() )
                {
                    pathIntermediate = pathShape;
                    return MS::kSuccess;
                }
            }
        }
        pathShape.pop();
    }
    return MS::kFailure;
}


MStatus cvWrapRebindCmd::createDriverSubsetMesh( MDagPath& pathDriver, MDagPath& pathDriverSubset )
{
    MStatus status;

    MFnMesh fnDriver( pathDriver );

    // Duplicate to create subset
    MStringArray duplicate;
    // Calling mesh,duplicate() gave jacked results.
    status = MGlobal::executeCommand( "duplicate -rr " + fnDriver.partialPathName(), duplicate );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    status = getDagPath( duplicate[0], pathDriverSubset );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    status = deleteIntermediateObjects( pathDriverSubset );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    // Invert selected faces so we can delete them
    int numFacesToDelete = fnDriver.numPolygons() - m_selectedFaces.length();
    if ( numFacesToDelete )
    {
        MIntArray facesToDelete( numFacesToDelete );
        int selectedFaceIndex = 0;
        int index = 0;
        for ( int i = 0; i < fnDriver.numPolygons(); i++ )
        {
            if ( i != m_selectedFaces[selectedFaceIndex] )
            {
                facesToDelete[index] = i;
                index++;
            }
            else
            {
                selectedFaceIndex++;
            }
        }

        MFnSingleIndexedComponent fnDeleteComp;
        MObject oFacesToDelete = fnDeleteComp.create( MFn::kMeshPolygonComponent, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = fnDeleteComp.addElements( facesToDelete );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        MSelectionList deleteList;
        status = deleteList.add( pathDriverSubset, oFacesToDelete );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = MGlobal::setActiveSelectionList( deleteList );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = MGlobal::executeCommand( "delete;" );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = getDagPath( duplicate[0], pathDriverSubset );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = pathDriverSubset.extendToShapeDirectlyBelow( 0 );
        CHECK_MSTATUS_AND_RETURN_IT( status );
    }
    return MS::kSuccess;
}


MStatus cvWrapRebindCmd::getMObject( MString& name, MObject& oNode )
{
    MStatus status;
    MSelectionList list;
    status = MGlobal::getSelectionListByName( name, list );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    status = list.getDependNode( 0, oNode );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    return MS::kSuccess;
}

MStatus cvWrapRebindCmd::getDagPath( MString& name, MDagPath& pathNode )
{
    MStatus status;
    MSelectionList list;
    status = MGlobal::getSelectionListByName( name, list );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    status = list.getDagPath( 0, pathNode );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    return MS::kSuccess;
}


MStatus cvWrapRebindCmd::deleteIntermediateObjects( MDagPath& pathNode )
{
    MStatus status;
    if ( pathNode.node().hasFn( MFn::kMesh ) )
    {
        status = pathNode.pop();
    }
    if ( !pathNode.node().hasFn( MFn::kTransform ) )
    {
        CHECK_MSTATUS_AND_RETURN_IT( MS::kFailure );
    }

    unsigned int childCount = pathNode.childCount();
    MFnDagNode fnTransform( pathNode );
    MString name = fnTransform.name();
    bool hasShape = false;
    for ( unsigned int i = 0; i < childCount; i++ )
    {
        MObject oChild = fnTransform.child( i, &status );
        if ( MFAIL( status ) )
        {
            continue;
        }
        status = pathNode.push( oChild );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        MFnDagNode fnNode( pathNode, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = pathNode.pop();
        CHECK_MSTATUS_AND_RETURN_IT( status );
        if ( fnNode.isIntermediateObject() )
        {
            status = MGlobal::executeCommand( "delete " + fnNode.partialPathName() );
        }
    }
    return MS::kSuccess;
}


MStatus cvWrapRebindCmd::undoIt()
{
    MStatus status;
    MObject oWrapNode, oCompData, oSampleIds, oDoubleData, oMatrixData, oNumericData;
    status = getMObject( m_wrapNode, oWrapNode );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MPlug plugSampleWeights( oWrapNode, cvWrap::aSampleWeights );
    MPlug plugSampleVerts( oWrapNode, cvWrap::aSampleComponents );
    MPlug plugSampleBindMatrix( oWrapNode, cvWrap::aBindMatrix );
    MPlug plugTriangleVerts( oWrapNode, cvWrap::aTriangleVerts );
    MPlug plugBarycentricWeights( oWrapNode, cvWrap::aBarycentricWeights );
    MPlug plugDirty( oWrapNode, cvWrap::aDirty );
    MFnSingleIndexedComponent fnComp;
    MFnComponentListData fnCompData;
    MFnMatrixData fnMatrixData;
    MFnNumericData fnNumericData;
    for ( unsigned int i = 0; i < m_selectedVerts.length(); i++ )
    {
        // Restore sample vert ids.
        if ( m_origSampleVerts.size() != 0 )
        {
            oSampleIds = fnComp.create( MFn::kMeshVertComponent, &status );
            CHECK_MSTATUS_AND_RETURN_IT( status );
            status = fnComp.addElements( m_origSampleVerts[i] );
            CHECK_MSTATUS_AND_RETURN_IT( status );
            oCompData = fnCompData.create( &status );
            CHECK_MSTATUS_AND_RETURN_IT( status );
            status = fnCompData.add( oSampleIds );
            CHECK_MSTATUS_AND_RETURN_IT( status );
            plugSampleVerts.elementByPhysicalIndex( m_selectedVerts[i], &status ).setMObject( oCompData );
            CHECK_MSTATUS_AND_RETURN_IT( status );
        }

        // Restore sample weights
        if ( m_origSampleWeights.size() != 0 )
        {
            MFnDoubleArrayData fnDoubleData;
            oDoubleData = fnDoubleData.create( m_origSampleWeights[i], &status );
            CHECK_MSTATUS_AND_RETURN_IT( status );
            plugSampleWeights.elementByPhysicalIndex( m_selectedVerts[i], &status ).setMObject( oDoubleData );
            CHECK_MSTATUS_AND_RETURN_IT( status );
        }

        // Restore bind matrix
        if ( m_origBindMatrices.length() != 0 )
        {
            oMatrixData = fnMatrixData.create( m_origBindMatrices[i], &status );
            CHECK_MSTATUS_AND_RETURN_IT( status );
            plugSampleBindMatrix.elementByPhysicalIndex( m_selectedVerts[i], &status ).setMObject( oMatrixData );
            CHECK_MSTATUS_AND_RETURN_IT( status );
        }

        // Restore triangle vertices
        if ( m_origTriangleVertices.size() != 0 )
        {
            oNumericData = fnNumericData.create( MFnNumericData::k3Int, &status );
            CHECK_MSTATUS_AND_RETURN_IT( status );
            status = fnNumericData.setData3Int( m_origTriangleVertices[i][0], m_origTriangleVertices[i][1], m_origTriangleVertices[i][2] );
            CHECK_MSTATUS_AND_RETURN_IT( status );
            plugTriangleVerts.elementByPhysicalIndex( m_selectedVerts[i], &status ).setMObject( oNumericData );
            CHECK_MSTATUS_AND_RETURN_IT( status );
        }

        // Store barycentric coordinates
        if ( m_origBarycentricCoordinates.size() != 0 )
        {
            oNumericData = fnNumericData.create( MFnNumericData::k3Float, &status );
            CHECK_MSTATUS_AND_RETURN_IT( status );
            status = fnNumericData.setData3Float( m_origBarycentricCoordinates[i][0], m_origBarycentricCoordinates[i][1], m_origBarycentricCoordinates[i][2] );
            CHECK_MSTATUS_AND_RETURN_IT( status );
            plugBarycentricWeights.elementByPhysicalIndex( m_selectedVerts[i], &status ).setMObject( oNumericData );
            CHECK_MSTATUS_AND_RETURN_IT( status );
        }
    }
    plugDirty.setBool( true );
    return MS::kSuccess;
}

