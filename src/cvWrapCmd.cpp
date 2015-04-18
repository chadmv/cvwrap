#include "cvWrapCmd.h"
#include "cvWrapDeformer.h"

#include <maya/MFnMatrixData.h>
#include <maya/MItSelectionList.h>

#define ZEROVECTOR(v) {v.x = 0.0; v.y = 0.0; v.z = 0.0;}
const char* CVWrapCmd::kName = "cvWrap";

/**
  Start a new progress bar.
  @param[in] title Status title.
  @param[in] count Progress bar maximum count.
*/
void StartProgress(const MString& title, unsigned int count) {
  if (MGlobal::mayaState() == MGlobal::kInteractive) {
    MString message = "progressBar -e -bp -ii true -st \"";
    message += title;
    message += "\" -max ";
    message += count;
    message += " $gMainProgressBar;";
    MGlobal::executeCommand(message);
  }
}


/**
  Increase the progress bar by the specified amount.
  @param[in] step Step amount.
*/
void StepProgress(int step) {
  if (MGlobal::mayaState() == MGlobal::kInteractive) {
    MString message = "progressBar -e -s ";
    message += step;
    message += " $gMainProgressBar;";
    MGlobal::executeCommand(message);
  }
}


/**
  Check if the progress has been cancelled.
  @return true if the progress has been cancelled.
*/
bool ProgressCancelled() {
  if (MGlobal::mayaState() == MGlobal::kInteractive) {
    int cmdResult = 0;
    MGlobal::executeCommand("progressBar -query -isCancelled $gMainProgressBar", cmdResult);
    return cmdResult != 0;
  }
  return false;
}


/**
  Ends any running progress bar.
*/
void EndProgress() {
  if (MGlobal::mayaState() == MGlobal::kInteractive) {
    MGlobal::executeCommand("progressBar -e -ep $gMainProgressBar;");
  }
}




CVWrapCmd::CVWrapCmd() {
}


MSyntax CVWrapCmd::newSyntax() {
  MSyntax syntax;
  syntax.addFlag("-n", "-name", MSyntax::kString);
  syntax.addFlag("-r", "-radius", MSyntax::kDouble);
  syntax.addFlag("-nbm", "-newBindMesh");
  syntax.addFlag("-ex", "-export", MSyntax::kString);
  syntax.addFlag("-im", "-import", MSyntax::kString);
  syntax.addFlag("-b", "-binding", MSyntax::kString);
  syntax.addFlag("-h", "-help");
  syntax.setObjectType(MSyntax::kSelectionList, 0, 255);
  syntax.useSelectionAsDefault(true);
  return syntax;
}


void* CVWrapCmd::creator() {																
  return new CVWrapCmd;										
}		

bool CVWrapCmd::isUndoable() const {
  return true;
}


MStatus CVWrapCmd::doIt(const MArgList& args) {
    MStatus status;
    m_radius = 0.1;
    m_name = "cvWrap#";
    m_import = false;
    m_export = false;
    m_useBinding = false;
    MArgDatabase argData( syntax(), args );
    m_newBindMesh = argData.isFlagSet( "-nbm" );
    m_displayHelp = argData.isFlagSet( "-h" );
    if ( argData.isFlagSet( "-r" ) )
    {
        m_radius = argData.flagArgumentDouble( "-r", 0, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        if ( m_radius <= 0.0 )
        {
            m_radius = 0.01;
        }
    }
    //m_radius *= m_radius;
    if ( argData.isFlagSet( "-n" ) )
    {
        m_name = argData.flagArgumentString( "-n", 0, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
    }
    if ( argData.isFlagSet( "-ex" ) )
    {
        m_export = true;
        m_filePath = argData.flagArgumentString( "-ex", 0, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
    }
    if ( argData.isFlagSet( "-im" ) )
    {
        m_import = true;
        m_filePath = argData.flagArgumentString( "-im", 0, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
    }
    if ( argData.isFlagSet( "-b" ) )
    {
        m_useBinding = true;
        m_filePath = argData.flagArgumentString( "-b", 0, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
    }

    if ( m_displayHelp )
    {
        displayHelp();
        return MS::kSuccess;
    }
    argData.getObjects( m_selectionList );

    MItSelectionList itSel( m_selectionList, MFn::kInvalid, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    if ( !m_import && !m_export )
    {
        status = getPathsToDriverAndDriven( m_pathDriver, m_pathDriven );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        MFnDagNode fnDriven( m_pathDriven );
        MString nameCmd = "-n \"" + m_name + "\" ";
        status = m_dgMod.commandToExecute( "deformer -type cvWrap " +
                  nameCmd + fnDriven.partialPathName() );
        CHECK_MSTATUS_AND_RETURN_IT( status );
    }
    else
    {
        status = m_selectionList.getDependNode( 0, m_oWrapNode );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        MFnDependencyNode fnNode( m_oWrapNode );
        if ( fnNode.typeId() != cvWrap::id )
        {
            MGlobal::displayError( "No wrap node specified." );
            return MS::kFailure;
        }
    }


    status = redoIt();
    CHECK_MSTATUS_AND_RETURN_IT( status );

    return MS::kSuccess;
}


MStatus CVWrapCmd::redoIt()
{
    MStatus status;


    if ( m_export || m_import )
    {
        if ( m_export )
        {
            status = exportBinding();
            CHECK_MSTATUS_AND_RETURN_IT( status );
        }
        else if ( m_import )
        {
            status = importBinding();
            CHECK_MSTATUS_AND_RETURN_IT( status );
        }
        return MS::kSuccess;
    }

    MFnDagNode fnDriver( m_pathDriver );
    MFnDagNode fnDriven( m_pathDriven );
    MFnMesh fnDriverMesh( m_pathDriver );
    MItGeometry itGeo( m_pathDriven, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MPointArray vertices;
    itGeo.allPositions( vertices );
    int numDrivenVerts = vertices.length();
    

    // Create the deformer
    status = m_dgMod.doIt();
    CHECK_MSTATUS_AND_RETURN_IT( status );
    // Reacquire the paths because on referenced geo, a new driven path is created.
    m_pathDriver = MDagPath();
    m_pathDriven = MDagPath();
    status = getPathsToDriverAndDriven( m_pathDriver, m_pathDriven );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MObject oDeformerNode;
    status = getLatestWrapNode( m_pathDriven, oDeformerNode );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    m_oWrapNode = oDeformerNode;
    MFnDependencyNode fnNode( oDeformerNode, &status );
    appendToResult( fnNode.name() );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    if ( m_useBinding )
    {
        status = importBinding();
        CHECK_MSTATUS_AND_RETURN_IT( status );
    }
    else
    {
        status = calculateBinding( oDeformerNode, m_pathDriver, m_pathDriven );
        CHECK_MSTATUS_AND_RETURN_IT( status );
    }

    

    // Connect the plugs
    MPlug plugDriverMesh = fnDriver.findPlug( "worldMesh", false, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    status = plugDriverMesh.selectAncestorLogicalIndex( 0, plugDriverMesh.attribute() );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    MPlug plugDriverGeo( oDeformerNode, cvWrap::aDriverGeo );
    MPlug plugBindMesh( oDeformerNode, cvWrap::aBindDriverGeo );

    // Create bind mesh so we can run rebind commands
    MDagPath pathBindMesh;
    MObject oDriver = m_pathDriver.node();
    status = getExistingBindMesh( oDriver, pathBindMesh );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    if ( m_newBindMesh || !pathBindMesh.isValid() )
    {
        MStringArray duplicate;
        // Calling mesh.duplicate() gives incorrect results.
        MGlobal::executeCommand( "duplicate -rr -n " + fnNode.name() + "Base " + fnDriver.partialPathName(), duplicate );
        status = getDagPath( duplicate[0], pathBindMesh );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = deleteIntermediateObjects( pathBindMesh );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        m_createdNodes.append( duplicate[0] );

        // Hide the duplicate
        MFnDagNode fnBindMesh( pathBindMesh, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        MPlug plug = fnBindMesh.findPlug( "visibility", &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = plug.setBool( false );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        status = pathBindMesh.extendToShapeDirectlyBelow( 0 );
        CHECK_MSTATUS_AND_RETURN_IT( status );

    }
    MFnDagNode fnBindMesh( pathBindMesh, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MPlug plugBindMessage = fnBindMesh.findPlug( "message", false, &status );   
    CHECK_MSTATUS_AND_RETURN_IT( status );

    // Store Nurbs surface attributes
    if ( m_pathDriven.node().apiType() == MFn::kNurbsSurface )
    {
        MFnNurbsSurface fnSurface( m_pathDriven, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        int cvsInU = fnSurface.numCVsInU();
        MPlug plugCVsInU( oDeformerNode, cvWrap::aCVsInU );
        plugCVsInU.setInt( cvsInU );

        int cvsInV = fnSurface.numCVsInV();
        MPlug plugCVsInV( oDeformerNode, cvWrap::aCVsInV );
        plugCVsInV.setInt( cvsInV );

        int formInU = fnSurface.formInU();
        MPlug plugFormInU( oDeformerNode, cvWrap::aFormInU );
        plugFormInU.setInt( formInU );

        int formInV = fnSurface.formInV();
        MPlug plugFormInV( oDeformerNode, cvWrap::aFormInV );
        plugFormInV.setInt( formInV );

        int degreeV = fnSurface.degreeV();
        MPlug plugDegreeV( oDeformerNode, cvWrap::aDegreeV );
        plugDegreeV.setInt( degreeV );
    }


    MDGModifier dgMod;
    dgMod.connect( plugDriverMesh, plugDriverGeo );
    dgMod.connect( plugBindMessage, plugBindMesh );
    status = dgMod.doIt();
    CHECK_MSTATUS_AND_RETURN_IT( status );


    return MS::kSuccess;
}


MStatus CVWrapCmd::getPathsToDriverAndDriven( MDagPath& pathDriver, MDagPath& pathDriven )
{
    MStatus status;
    MItSelectionList itSel( m_selectionList, MFn::kInvalid, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    MObject oComponent;
    MDagPath tempDag;
    for ( ; !itSel.isDone(); itSel.next() )
    {
        itSel.getDagPath( tempDag, oComponent );
        unsigned int childCount = tempDag.childCount();
        
        for ( unsigned int j = 0; j < childCount; j++ )
        {
            MObject oChild = tempDag.child( j, &status );
            CHECK_MSTATUS_AND_RETURN_IT( status );
            tempDag.push( oChild );

            MFnDagNode fnDag( tempDag );
            if ( fnDag.isIntermediateObject() )
            {
                tempDag.pop();
                continue;
            }

            if ( tempDag.node().apiType() == MFn::kMesh || tempDag.node().apiType() == MFn::kNurbsSurface || tempDag.node().apiType() == MFn::kNurbsCurve )
            {
                if ( !pathDriven.isValid() )
                {
                    pathDriven = tempDag;
                    break;
                }
                else
                {
                    pathDriver = tempDag;
                    break;
                }
            }
        }
    }
    return MS::kSuccess;
}


MStatus CVWrapCmd::getLatestWrapNode( MDagPath& pathDriven, MObject& oDeformerNode )
{
    MStatus status;
    MObject oDriven = pathDriven.node();
    MItDependencyGraph itDG( 
        oDriven, 
        MFn::kGeometryFilt, 
        MItDependencyGraph::kUpstream, 
        MItDependencyGraph::kDepthFirst,
        MItDependencyGraph::kNodeLevel, 
        &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    for ( ; !itDG.isDone(); itDG.next() )
    {
        oDeformerNode = itDG.currentItem();
        MFnDependencyNode fnNode( oDeformerNode, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        if ( fnNode.typeId() == cvWrap::id )
        {
            return MS::kSuccess;
        }
    }
    return MS::kFailure;
}


MStatus CVWrapCmd::calculateBinding( MObject& oWrapNode, MDagPath& pathDriver, MDagPath& pathDriven )
{
    MStatus status;
    MMeshIntersector intersector;
    MObject oDriver = pathDriver.node();
    MMatrix driverMatrix = pathDriver.inclusiveMatrix();
    status = intersector.create( oDriver, driverMatrix );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    MFnMesh fnDriverMesh( pathDriver, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    std::vector<std::set<int> > adjacency;
    status = getAdjacency( pathDriver, adjacency );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    MPointArray points;
    MFloatVectorArray normals;
    fnDriverMesh.getPoints( points, MSpace::kWorld );
    fnDriverMesh.getVertexNormals( false, normals, MSpace::kWorld );

    MPoint pt, closestPoint;
    MPointOnMesh pointOnMesh;
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
    MPlug plugBindInfoElement;
    MFnNumericData fnNumericData;
    MFnMatrixData fnMatrixData;
    MMatrix matrix;
    MItGeometry itGeo( pathDriven, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    int triangleVertices[3];
    float a, b, c;
    double totalWeight;


    bool progressbar = MGlobal::mayaState() == MGlobal::kInteractive;
    // Start progress bar
    if (progressbar) {
      StartProgress("Binding wrap...", itGeo.count());
    }

    MPointArray inputPoints;
    itGeo.allPositions( inputPoints, MSpace::kWorld );

    int ii = 0;
    std::map<int, double>::iterator itDistance;
    for ( itGeo.reset(); !itGeo.isDone(); itGeo.next(), ii++ )
    {

        // Get closest point and faceId
        pt = inputPoints[ii];
        status = intersector.getClosestPoint( pt, pointOnMesh );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        closestPoint = MPoint( pointOnMesh.getPoint() ) * driverMatrix;
        faceId = pointOnMesh.faceIndex();
        triangleId = pointOnMesh.triangleIndex();

        // Get vertices of closest triangle
        vertexList.clear();
        status = fnDriverMesh.getPolygonTriangleVertices( faceId, triangleId, triangleVertices );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Get barycentric coordinates of closestPoint
        getBarycentricCoordinates( closestPoint,
                points[triangleVertices[0]],
                points[triangleVertices[1]],
                points[triangleVertices[2]],
                a, b, c );

        // Get vertices of closest face
        vertexList.clear();
        status = fnDriverMesh.getPolygonVertices( faceId, vertexList );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Crawl surface to get sample weights
        totalWeight = 0.0;
        std::map<int, double> distances;
        for ( unsigned int i = 0; i < vertexList.length(); ++i )
        {
            status = crawlSurface( vertexList[i], points, distances,
                    0.0, closestPoint, m_radius, adjacency );
            CHECK_MSTATUS_AND_RETURN_IT( status );
        }

        calculateSampleWeights( distances, sampleIds, weights, normalizedWeights,
                closestPoint, points, totalWeight );

        if ( distances.size() < vertexList.length() || totalWeight == 0.0 )
        {
            distances.clear();
            for ( unsigned int i = 0; i < vertexList.length(); ++i )
            {
                distances[vertexList[i]] = closestPoint.distanceTo( points[vertexList[i]] );
                //distances[vertexList[i]] = distanceSquared( closestPoint, points[vertexList[i]] );
            }
            calculateSampleWeights( distances, sampleIds, weights, normalizedWeights,
                    closestPoint, points, totalWeight );
            if ( totalWeight == 0.0 )
            {
                // Just take equal weight from each point
                sampleIds.setLength( vertexList.length() );
                weights.setLength( vertexList.length() );
                normalizedWeights.setLength( vertexList.length() );
                double w = 1.0 / vertexList.length();
                for ( unsigned int i = 0; i < vertexList.length(); ++i )
                {
                    sampleIds[i] = vertexList[i];
                    weights[i] = w;
                    normalizedWeights[i] = w;
                }
            }
        }

        // Calculate origin
        origin = closestPoint;

        // Calculate normal and up
        normal = MVector::zero;
        up = MVector::zero;
        for ( unsigned int i = 0; i < weights.length(); i++ )
        {
            normal += MVector( normals[sampleIds[i]] ) * normalizedWeights[i];
            up += (points[sampleIds[i]] - origin) * normalizedWeights[i];
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
                up -= (points[sortedWeights[i].index] - origin) * sortedWeights[i].normalizedWeight;
            }
        }
        up.normalize();
        normal.normalize();


        // Store sample vert ids.
        oComp = fnComp.create( MFn::kMeshVertComponent, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = fnComp.addElements( sampleIds );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        oCompData = fnCompData.create();
        fnCompData.add( oComp );
        plugSampleVerts.elementByLogicalIndex( itGeo.index(), &status ).setMObject( oCompData );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Store sample weights
        MFnDoubleArrayData fnDoubleData;
        oDoubleData = fnDoubleData.create( weights, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        plugSampleWeights.elementByLogicalIndex( itGeo.index(), &status ).setMObject( oDoubleData );
        CHECK_MSTATUS_AND_RETURN_IT( status );


        // Store bind matrix
        cvWrap::createMatrix( matrix, origin, normal, up );
        matrix = matrix.inverse();
        oMatrixData = fnMatrixData.create( matrix, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        plugSampleBindMatrix.elementByLogicalIndex( itGeo.index(), &status ).setMObject( oMatrixData );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Store triangle vertices
        oNumericData = fnNumericData.create( MFnNumericData::k3Int, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = fnNumericData.setData3Int( triangleVertices[0], triangleVertices[1], triangleVertices[2] );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        plugTriangleVerts.elementByLogicalIndex( itGeo.index(), &status ).setMObject( oNumericData );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Store barycentric coordinates
        oNumericData = fnNumericData.create( MFnNumericData::k3Float, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = fnNumericData.setData3Float( a, b, c );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        plugBarycentricWeights.elementByLogicalIndex( itGeo.index(), &status ).setMObject( oNumericData );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Increment the progress bar on every 250 verts because updating the UI is slow.
        if (itGeo.index() % 250 == 0 && itGeo.index() != 0 ) {
          StepProgress(250);
          if (ProgressCancelled()) {
            break;
          }
        }
    }

    EndProgress();
    return MS::kSuccess;
}


void CVWrapCmd::getBarycentricCoordinates( MPoint& P, MPoint& A, MPoint& B, MPoint& C, float& a, float& b, float& c )
{
    // Compute the normal of the triangle
    MVector N = (B - A) ^ (C - A);
    MVector unitN = N.normal();

    // Compute twice area of triangle ABC
    double areaABC = unitN * N;

    if ( areaABC == 0.0f )
    {
        a = 0.33f;
        b = 0.33f;
        c = 0.33f;
        return;
    }

    // Compute a
    double areaPBC = unitN * ((B - P) ^ (C - P));
    a = (float)(areaPBC / areaABC);

    // Compute b
    double areaPCA = unitN * ((C - P) ^ (A - P));
    b = (float)(areaPCA / areaABC);

    // Compute c
    c = 1.0f - a - b;
}


MStatus CVWrapCmd::getAdjacency( MDagPath& pathMesh, std::vector<std::set<int> >& adjacency )
{
    MStatus status;
    // Get mesh adjacency
    MIntArray faces, vertices;
    MObject oNull = MObject::kNullObj;
    MItMeshVertex itVert( pathMesh, oNull, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MFnMesh fnMesh( pathMesh, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    adjacency.resize( itVert.count() );
    for ( ; !itVert.isDone(); itVert.next() )
    {
        faces.clear();
        status = itVert.getConnectedFaces( faces );
        adjacency[itVert.index()].clear();
        for ( unsigned int j = 0; j < faces.length(); ++j )
        {
            vertices.clear();
            fnMesh.getPolygonVertices( faces[j], vertices );
            for ( unsigned int k = 0; k < vertices.length(); ++k )
            {
                if ( vertices[k] != itVert.index() )
                {
                    adjacency[itVert.index()].insert( vertices[k] );
                }
            }
        }
    }
    return MS::kSuccess;
}


/**
    @brief Recursive function that crawls the surface to find all the points
    within the projected brush radius.

    @param[in] geom_index Geometry index.
    @param[in] vertex_index Current vertex index to crawl.
    @param[out] distances Storage for the distances of the crawled points.
    @param[in] source_distance Previous calculated distance.
    @param[in] source_point Previous crawled point.

    @return MStatus
 */
MStatus CVWrapCmd::crawlSurface( int vertexIndex,
        MPointArray& points, std::map<int, double>& distances, double sourceDistance,
        MPoint& sourcePoint, double maxDistance, std::vector<std::set<int> >& adjacency )
{
    MStatus status;
    if ( distances.size() > 50 )
    {
        return MS::kSuccess;
    }
    if ( sourceDistance >= maxDistance )
    {
        return MS::kSuccess;
    }

    MPoint& pt = points[vertexIndex];
    sourceDistance += sourcePoint.distanceTo( pt );
    //sourceDistance += distanceSquared( sourcePoint, pt );
    if ( sourceDistance >= maxDistance )
    {
        return MS::kSuccess;
    }
    if ( sourceDistance <= distances[vertexIndex] || distances[vertexIndex] < 0.0001 )
    {
        distances[vertexIndex] = sourceDistance;
    }
    else
    {
        // A smaller distance is already stored so we don't want to crawl
        // from this vertex any further.
        return MS::kSuccess;
    }

    // Crawl the connected vertices
    std::set<int>::iterator itAdjacency;
    for ( itAdjacency = adjacency[vertexIndex].begin();
            itAdjacency != adjacency[vertexIndex].end();
            itAdjacency++ )
    {
        status = crawlSurface( *itAdjacency, points, 
                distances, sourceDistance, pt, maxDistance, adjacency );
        CHECK_MSTATUS_AND_RETURN_IT( status );
    }

    return MS::kSuccess;
}


double CVWrapCmd::distanceSquared( const MPoint& p1, const MPoint& p2 )
{
    double xx = p2.x - p1.x;
    double yy = p2.y - p1.y;
    double zz = p2.z - p1.z;
    return (xx*xx) + (yy*yy) + (zz*zz);
}


void CVWrapCmd::calculateSampleWeights( std::map<int, double>& distances,
        MIntArray& sampleIds,
        MDoubleArray& weights,
        MDoubleArray& normalizedWeights,
        MPoint& closestPoint,
        MPointArray& points,
        double& totalWeight )
{
    std::map<int, double>::iterator itDistance;
    unsigned int length = (unsigned int)distances.size();
    sampleIds.setLength(length);
    weights.setLength(length);
    normalizedWeights.setLength(length);
    int ii = 0;
    totalWeight = 0.0;
    double maxDistance = 0.0;
    for ( itDistance = distances.begin();
            itDistance != distances.end();
            itDistance++, ii++ )
    {
        sampleIds[ii] = itDistance->first;
        weights[ii] = closestPoint.distanceTo( points[sampleIds[ii]] );
        //weights[ii] = distanceSquared( closestPoint, points[sampleIds[ii]] );
        if ( weights[ii] > maxDistance )
        {
            maxDistance = weights[ii];
        }
    }
    ii = 0;
    for ( itDistance = distances.begin();
            itDistance != distances.end();
            itDistance++, ii++ )
    {
        weights[ii] = 1.0 - (itDistance->second / maxDistance);
        totalWeight += weights[ii];
    }
    if ( totalWeight == 0.0 )
    {
        return;
    }
    for ( unsigned int i = 0; i < weights.length(); i++ )
    {
        normalizedWeights[i] = weights[i] / totalWeight;
    }
}


/**
    @brief Gets the MDagPath of any existing bind wrap mesh so we don't have to duplicate it for each new wrap.

    @param[in] oDriver Driver mesh MObject
    @param[out] pathBindMesh Storage for path to an existing bind mesh
*/
MStatus CVWrapCmd::getExistingBindMesh( MObject& oDriver, MDagPath &pathBindMesh )
{
    MStatus status;
    MFnDependencyNode fnDriver( oDriver, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    MPlug plugOutGeom = fnDriver.findPlug( "outMesh", false, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MPlugArray geomPlugs;
    plugOutGeom.connectedTo( geomPlugs, false, true );
    for ( unsigned int i = 0; i < geomPlugs.length(); i++ )
    {
        MObject oThisNode = geomPlugs[i].node();
        MFnDependencyNode fnNode( oThisNode );
        if ( fnNode.typeId() == cvWrap::id )
        {
            // Get bind wrap mesh from wrap node
            MPlug plugBindWrapMesh = fnNode.findPlug( "bindMesh", false, &status );
            CHECK_MSTATUS_AND_RETURN_IT( status );
            MPlugArray bindPlugs;
            plugBindWrapMesh.connectedTo( bindPlugs, true, false );
            if ( bindPlugs.length() > 0 )
            {
                MObject oBindMesh = bindPlugs[0].node();
                MFnDagNode fnBindDag( oBindMesh, &status );
                CHECK_MSTATUS_AND_RETURN_IT( status );
                status = fnBindDag.getPath( pathBindMesh );
                CHECK_MSTATUS_AND_RETURN_IT( status );
            }
            break;
        }
    }
    return MS::kSuccess;
}



MStatus CVWrapCmd::getDagPath( MString& name, MDagPath& pathNode )
{
    MStatus status;
    MSelectionList list;
    status = MGlobal::getSelectionListByName( name, list );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    status = list.getDagPath( 0, pathNode );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    return MS::kSuccess;
}


MStatus CVWrapCmd::deleteIntermediateObjects( MDagPath& pathNode )
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


MStatus CVWrapCmd::undoIt()
{
    MStatus status;
    status = m_dgMod.undoIt();
    CHECK_MSTATUS_AND_RETURN_IT( status );

    for ( unsigned int i = 0; i < m_createdNodes.length(); i++ )
    {
        MGlobal::executeCommand( "delete " + m_createdNodes[i] );
    }
    m_createdNodes.clear();

    return MS::kSuccess;
}


MStatus CVWrapCmd::exportBinding()
{
    MStatus status;
    MFnDependencyNode fnWrapNode( m_oWrapNode, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    if ( fnWrapNode.typeId() != cvWrap::id )
    {
        CHECK_MSTATUS_AND_RETURN_IT( MS::kFailure );
    }

    std::ofstream out( m_filePath.asChar(), ios::binary );
    if ( !out.is_open() )
    {
        MGlobal::displayError( "Unable to open file for writing." );
        CHECK_MSTATUS_AND_RETURN_IT( MS::kFailure );
    }

    MFnSingleIndexedComponent fnComp;
    MFnComponentListData fnCompData;
    MFnNumericData fnNumericData;
    MFnMatrixData fnMatrixData;
    MMatrix matrix;
    MObject oTemp;
    MFloatArray tempCoords( 3 );
    MIntArray triangleVerts( 3 );
    MDoubleArray tempWeights;
    MIntArray sampleIds;
    MIntArray tempIds;
    MObject oComp;

    MPlug plugSampleWeights( m_oWrapNode, cvWrap::aSampleWeights );
    MPlug plugSampleVerts( m_oWrapNode, cvWrap::aSampleComponents );
    MPlug plugSampleBindMatrix( m_oWrapNode, cvWrap::aBindMatrix );
    MPlug plugTriangleVerts( m_oWrapNode, cvWrap::aTriangleVerts );
    MPlug plugBarycentricWeights( m_oWrapNode, cvWrap::aBarycentricWeights );

    float version = 1.0f;
    out.write( (char *)&version, sizeof( float ) );

    unsigned int numVerts = plugSampleWeights.numElements();
    out.write( (char *)(&numVerts), sizeof( int ) );

    // Start progress bar
    StartProgress("Exporting wrap...", numVerts);
    for (unsigned int i = 0; i < numVerts; ++i) {
        // Export sample vertex ids
        oTemp = plugSampleVerts.elementByPhysicalIndex( i, &status ).asMObject();
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = fnCompData.setObject( oTemp );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        sampleIds.clear();
        for ( unsigned int j = 0; j < fnCompData.length(); j++ )
        {
            oComp = fnCompData[j];
            if ( oComp.hasFn( MFn::kSingleIndexedComponent ) )
            {
                status = fnComp.setObject( oComp );
                CHECK_MSTATUS_AND_RETURN_IT( status );
                fnComp.getElements( tempIds );
                for( unsigned int k = 0; k < tempIds.length(); k++ )
                {
                    sampleIds.append( tempIds[k] );
                }
            }
        }
        writeAttribute( out, sampleIds );

        // Export sample weights
        oTemp = plugSampleWeights.elementByPhysicalIndex( i, &status ).asMObject();
        CHECK_MSTATUS_AND_RETURN_IT( status );
        MFnDoubleArrayData fnDoubleData( oTemp );
        tempWeights = fnDoubleData.array();
        writeAttribute( out, tempWeights );

        // Export bind matrix
        oTemp = plugSampleBindMatrix.elementByPhysicalIndex( i, &status ).asMObject();
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = fnMatrixData.setObject( oTemp );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        writeAttribute( out, fnMatrixData.matrix() );

        // Export triangle vertices
        oTemp = plugTriangleVerts.elementByPhysicalIndex( i, &status ).asMObject();
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = fnNumericData.setObject( oTemp );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        fnNumericData.getData3Int( triangleVerts[0], triangleVerts[1], triangleVerts[2] );
        writeAttribute( out, triangleVerts );

        // Export barycentric coordinates
        oTemp = plugBarycentricWeights.elementByPhysicalIndex( i, &status ).asMObject();
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = fnNumericData.setObject( oTemp );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        fnNumericData.getData3Float( tempCoords[0], tempCoords[1], tempCoords[2] );
        writeAttribute( out, tempCoords );

        // Increment the progress bar
        if (i % 250 == 0 && i != 0) {
          StepProgress(250);
          if(ProgressCancelled()) {
            break;
          }
        }
    }
    EndProgress();
    out.close();

    MGlobal::displayInfo( "Wrap binding exported." );

    return status;
}


MStatus CVWrapCmd::writeAttribute( std::ofstream &out, const MIntArray &attribute )
{
    unsigned int length = attribute.length();
    out.write( (char *)(&length), sizeof( int ) );
    if ( length > 0 )
    {
        int * pAttr = new int[length];
        attribute.get( pAttr );
        out.write( (char *)pAttr, length * sizeof( int ) );
        delete [] pAttr;
    }
    return MS::kSuccess;
}


MStatus CVWrapCmd::writeAttribute( std::ofstream &out, const MDoubleArray &attribute )
{
    unsigned int length = attribute.length();
    out.write( (char *)(&length), sizeof( int ) );
    if ( length > 0 )
    {
        double * pAttr = new double[length];
        attribute.get( pAttr );
        out.write( (char *)pAttr, length * sizeof( double ) );
        delete [] pAttr;
    }
    return MS::kSuccess;
}


MStatus CVWrapCmd::writeAttribute( std::ofstream &out, const MFloatArray &attribute )
{
    unsigned int length = attribute.length();
    out.write( (char *)(&length), sizeof( int ) );
    if ( length > 0 )
    {
        float * pAttr = new float[length];
        attribute.get( pAttr );
        out.write( (char *)pAttr, length * sizeof( float ) );
        delete [] pAttr;
    }
    return MS::kSuccess;
}


MStatus CVWrapCmd::writeAttribute( std::ofstream &out, const MMatrix &attribute )
{
    unsigned int length = 16;
    out.write( (char *)(&length), sizeof( int ) );
    double pAttr[16];
    for ( int i = 0; i < 4; i++ )
    {
        for ( int j = 0; j < 4; j++ )
        {
            pAttr[(i * 4) + j] = attribute[i][j];
        }
    }
    out.write( (char *)pAttr, length * sizeof( double ) );
    return MS::kSuccess;
}


MStatus CVWrapCmd::importBinding()
{
    MStatus status;

    std::ifstream in( m_filePath.asChar(), ios::binary );
    if ( !in.is_open() )
    {
        MGlobal::displayInfo( "Unable to open file for importing." );
        CHECK_MSTATUS_AND_RETURN_IT( MS::kFailure );
    }

    MFnSingleIndexedComponent fnComp;
    MFnComponentListData fnCompData;
    MFnNumericData fnNumericData;
    MFnMatrixData fnMatrixData;
    MMatrix matrix;
    MObject oTemp, oData;
    MPlug plugTemp;
    MFloatArray tempCoords( 3 );
    MIntArray triangleVerts( 3 );
    MDoubleArray tempWeights;
    MIntArray sampleIds;
    MFnDoubleArrayData fnDoubleData;

    MPlug plugSampleWeights( m_oWrapNode, cvWrap::aSampleWeights );
    MPlug plugSampleVerts( m_oWrapNode, cvWrap::aSampleComponents );
    MPlug plugSampleBindMatrix( m_oWrapNode, cvWrap::aBindMatrix );
    MPlug plugTriangleVerts( m_oWrapNode, cvWrap::aTriangleVerts );
    MPlug plugBarycentricWeights( m_oWrapNode, cvWrap::aBarycentricWeights );

    float version;
    in.read( (char *)(&version), sizeof( float ) );

    int numVerts;
    in.read( (char *)(&numVerts), sizeof( int ) );
    //if ( numVerts != numDrivenVerts )
    //{
    //    char buffer[256];
    //    sprintf( buffer, "Mesh has %d vertices. Binding file containts %d vertices", numDrivenVerts, numVerts );
    //    MGlobal::displayError( buffer );
    //    in.close();
    //    return MS::kFailure;
    //}
    StartProgress("Importing wrap...", numVerts);
    for ( int i = 0; i < numVerts; i++ )
    {
        // Import sample vertex ids
        readAttribute( in, sampleIds );
        oTemp = fnComp.create( MFn::kMeshVertComponent, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = fnComp.addElements( sampleIds );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        oData = fnCompData.create( &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = fnCompData.add( oTemp );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        plugSampleVerts.elementByLogicalIndex( i, &status ).setMObject( oData );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Import sample weights
        readAttribute( in, tempWeights );
        MFnDoubleArrayData fnDoubleData;
        oData = fnDoubleData.create( tempWeights, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        plugSampleWeights.elementByLogicalIndex( i, &status ).setMObject( oData );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Import bind matrix
        readAttribute( in, matrix );
        oData = fnMatrixData.create( matrix );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        plugSampleBindMatrix.elementByLogicalIndex( i, &status ).setMObject( oData );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Import triangle vertices
        readAttribute( in, triangleVerts );
        oData = fnNumericData.create( MFnNumericData::k3Int, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = fnNumericData.setData3Int( triangleVerts[0], triangleVerts[1], triangleVerts[2] );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        plugTriangleVerts.elementByLogicalIndex( i, &status ).setMObject( oData );
        CHECK_MSTATUS_AND_RETURN_IT( status );

        // Import barycentric coordinates
        readAttribute( in, tempCoords );
        oData = fnNumericData.create( MFnNumericData::k3Float, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        fnNumericData.setData3Float( tempCoords[0], tempCoords[1], tempCoords[2] );
        plugBarycentricWeights.elementByLogicalIndex( i, &status ).setMObject( oData );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        // Increment the progress bar
        if (i % 250 == 0 && i != 0) {
          StepProgress(250);
          if (ProgressCancelled()) {
            break;
          }
        }
    }
    EndProgress();
    in.close();
    MGlobal::displayInfo( "Wrap binding imported." );

    return MS::kSuccess;
}


MStatus CVWrapCmd::readAttribute( std::ifstream &in, MIntArray &attribute )
{
    attribute.clear();
    int numValues;
    in.read( (char *)(&numValues), sizeof( int ) );
    if ( numValues > 0 )
    {
        attribute.setLength( numValues );
        int * pValues = new int[numValues];
        in.read( (char *)pValues, numValues * sizeof( int ) );
        for ( int i = 0; i < numValues; i++ )
        {
            attribute[i] = pValues[i];
        }
        delete [] pValues;
    }

    return MS::kSuccess;
}


MStatus CVWrapCmd::readAttribute( std::ifstream &in, MDoubleArray &attribute )
{
    attribute.clear();
    int numValues;
    in.read( (char *)(&numValues), sizeof( int ) );
    if ( numValues > 0 )
    {
        attribute.setLength( numValues );
        double * pValues = new double[numValues];
        in.read( (char *)pValues, numValues * sizeof( double ) );
        for ( int i = 0; i < numValues; i++ )
        {
            attribute[i] = pValues[i];
        }
        delete [] pValues;
    }
    return MS::kSuccess;
}


MStatus CVWrapCmd::readAttribute( std::ifstream &in, MFloatArray &attribute )
{
    attribute.clear();
    int numValues;
    in.read( (char *)(&numValues), sizeof( int ) );
    if ( numValues > 0 )
    {
        attribute.setLength( numValues );
        float * pValues = new float[numValues];
        in.read( (char *)pValues, numValues * sizeof( float ) );
        for ( int i = 0; i < numValues; i++ )
        {
            attribute[i] = pValues[i];
        }
        delete [] pValues;
    }
    return MS::kSuccess;
}


MStatus CVWrapCmd::readAttribute( std::ifstream &in, MMatrix &matrix )
{
    int numValues;
    in.read( (char *)(&numValues), sizeof( int ) );
    if ( numValues > 0 )
    {
        double * pValues = new double[numValues];
        in.read( (char *)pValues, numValues * sizeof( double ) );
        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                matrix[i][j] = pValues[(i * 4) + j];
            }
        }
        delete [] pValues;
    }
    return MS::kSuccess;
}


void CVWrapCmd::displayHelp()
{
   MString help;
   help += "Flags:\n"; 
   help += "-name (-n):          String     Name of the wrap node to create.\n"; 
   help += "-radius (-r):        Double     Sample radius.  Default is 0.1.  The greater the radius,\n"; 
   help += "                                the smoother the deformation but slower performance.\n";
   help += "-newBindMesh (-nbm)  N/A        Creates a new bind mesh, otherwise the existing bind mesh will be used.\n";
   help += "-export (-ex):       String     Path to a file to export the binding to.\n"; 
   help += "-import (-im):       String     Path to a file to import the binding from.\n"; 
   help += "-binding (-b):       String     Path to a file to import the binding from on creation.\n"; 
   help += "-help (-h)           N/A        Display this text.\n";
   MGlobal::displayInfo( help );
}
