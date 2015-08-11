#include "bindingio.h"
#include "cvWrapDeformer.h"

#include <maya/MGlobal.h>
#include <maya/MObjectArray.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MFnMatrixData.h>
#include <maya/MFnWeightGeometryFilter.h>

const float BindingIO::kWrapFileVersion = 1.0f;

template <>
void WriteAttribute<double, MMatrix>(std::ofstream &out, const MMatrix& attribute) {
  double values[16];
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      values[i*4 + j] = attribute[i][j];
    }
  }
  out.write((char *)values, 16 * sizeof(double));
}

template <>
void ReadAttribute<double, MMatrix>(std::ifstream &in, MMatrix &matrix) {
  double values[16];
  in.read((char *)values, 16 * sizeof(double));
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      matrix[i][j] = values[(i * 4) + j];
    }
  }
}

MStatus BindingIO::ExportBinding(std::ofstream& out, MObject& oWrapNode) {
  MStatus status;
  MFnWeightGeometryFilter fnWrapNode(oWrapNode, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  if (fnWrapNode.typeId() != CVWrap::id) {
    MGlobal::displayError(fnWrapNode.name() + " is not a cvWrap node.");
    CHECK_MSTATUS_AND_RETURN_IT(MS::kFailure);
  }

  out.write((char *)&kWrapFileVersion, sizeof(float));

  MPlug plugBindData = fnWrapNode.findPlug(CVWrap::aBindData, false, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  // Get the input geometry so we can get the geometry indices
  MObjectArray outputGeometry;
  status = fnWrapNode.getOutputGeometry(outputGeometry);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  // Write the number of geometry
  unsigned int geometryCount = outputGeometry.length();
  out.write((char *)(&geometryCount), sizeof(geometryCount));

  MIntArray triangleVerts(3);  /**< Storage for the triangle vertex ids. */
  MFloatArray baryCoords(3);  /**< Storage for the barycentric weights. */
  for (unsigned int i = 0; i < outputGeometry.length(); ++i) {
    unsigned int geomIndex = fnWrapNode.indexForOutputShape(outputGeometry[i], &status);
    // Get the plugs to the binding attributes for this geometry
    MPlug plugBind = plugBindData.elementByLogicalIndex(geomIndex, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MPlug plugSampleWeights = plugBind.child(CVWrap::aSampleWeights, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MPlug plugSampleVerts = plugBind.child(CVWrap::aSampleComponents, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MPlug plugSampleBindMatrix = plugBind.child(CVWrap::aBindMatrix, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MPlug plugTriangleVerts = plugBind.child(CVWrap::aTriangleVerts, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MPlug plugBarycentricWeights = plugBind.child(CVWrap::aBarycentricWeights, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);

    unsigned int numElements = plugSampleWeights.numElements();
    out.write((char *)(&numElements), sizeof(numElements));
  
    for (unsigned int j = 0; j < numElements; ++j) {
      // Write the logical index
      MPlug plugSampleVertElement = plugSampleVerts.elementByPhysicalIndex(j, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      unsigned int logicalIndex = plugSampleVertElement.logicalIndex();
      out.write((char *)(&logicalIndex), sizeof(logicalIndex));

      // Export sample vertex ids
      MObject oSampleIds = plugSampleVertElement.asMObject();
      CHECK_MSTATUS_AND_RETURN_IT(status);
      MFnIntArrayData fnIntData(oSampleIds, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      MIntArray sampleIds = fnIntData.array();
      WriteAttribute<int, MIntArray>(out, sampleIds);

      // Export sample weights
      MObject oWeightData = plugSampleWeights.elementByPhysicalIndex(j, &status).asMObject();
      CHECK_MSTATUS_AND_RETURN_IT(status);
      MFnDoubleArrayData fnDoubleData(oWeightData);
      MDoubleArray weights = fnDoubleData.array();
      WriteAttribute<double, MDoubleArray>(out, weights);

      // Export bind matrix
      MObject oBindMatrix = plugSampleBindMatrix.elementByPhysicalIndex(j, &status).asMObject();
      CHECK_MSTATUS_AND_RETURN_IT(status);
      MFnMatrixData fnMatrixData(oBindMatrix, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      WriteAttribute<double, MMatrix>(out, fnMatrixData.matrix());

      // Export triangle vertices
      MObject oTriangleVerts = plugTriangleVerts.elementByPhysicalIndex(j, &status).asMObject();
      CHECK_MSTATUS_AND_RETURN_IT(status);
      MFnNumericData fnNumericData(oTriangleVerts, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      fnNumericData.getData3Int(triangleVerts[0], triangleVerts[1], triangleVerts[2]);
      WriteAttribute<int, MIntArray>(out, triangleVerts);

      // Export the barycentric weights
      MObject oBaryWeights = plugBarycentricWeights.elementByPhysicalIndex(j, &status).asMObject();
      CHECK_MSTATUS_AND_RETURN_IT(status);
      MFnNumericData fnBaryData(oBaryWeights, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      fnBaryData.getData3Float(baryCoords[0], baryCoords[1], baryCoords[2]);
      WriteAttribute<float, MFloatArray>(out, baryCoords);
    }
  }

  MGlobal::displayInfo("Wrap binding exported.");

  return status;
}



MStatus BindingIO::ImportBinding(std::ifstream& in, MObject& oWrapNode) {
  MStatus status;

  MFnWeightGeometryFilter fnWrapNode(oWrapNode, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  MPlug plugBindData = fnWrapNode.findPlug(CVWrap::aBindData, false, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  float version;
  in.read((char *)(&version), sizeof(float));

  unsigned int geometryCount = 0;
  in.read((char *)(&geometryCount), sizeof(geometryCount));

  MFnMatrixData fnMatrixData;
  MFnIntArrayData fnIntData;
  MFnDoubleArrayData fnDoubleData;
  MFnNumericData fnNumericData;
  // We are assuming that the geometryIndices are compact and continuous.  It is possible
  // that the indices could be sparse, but we will ignore that corner case.
  for (unsigned int geomIndex = 0; geomIndex < geometryCount; ++geomIndex) {
    // Get the plugs to the binding attributes for this geometry
    MPlug plugBind = plugBindData.elementByLogicalIndex(geomIndex, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MPlug plugSampleWeights = plugBind.child(CVWrap::aSampleWeights, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MPlug plugSampleVerts = plugBind.child(CVWrap::aSampleComponents, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MPlug plugSampleBindMatrix = plugBind.child(CVWrap::aBindMatrix, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MPlug plugTriangleVerts = plugBind.child(CVWrap::aTriangleVerts, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    MPlug plugBarycentricWeights = plugBind.child(CVWrap::aBarycentricWeights, &status);
    CHECK_MSTATUS_AND_RETURN_IT(status);

    unsigned int numElements = plugSampleWeights.numElements();
    in.read((char *)(&numElements), sizeof(numElements));
    for (unsigned int i = 0; i < numElements; ++i) {
      unsigned int logicalIndex = 0;
      in.read((char *)(&logicalIndex), sizeof(logicalIndex));

      // Sample vert ids.
      MIntArray sampleIds;
      ReadAttribute<int, MIntArray>(in, sampleIds);
      MObject oIntData = fnIntData.create(sampleIds, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      plugSampleVerts.elementByLogicalIndex(logicalIndex, &status).setMObject(oIntData);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      // Sample weights
      MDoubleArray weights;
      ReadAttribute<double, MDoubleArray>(in, weights);
      MObject oDoubleData = fnDoubleData.create(weights, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      plugSampleWeights.elementByLogicalIndex(logicalIndex, &status).setMObject(oDoubleData);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      // Bind matrix
      MMatrix bindMatrix;
      ReadAttribute<double, MMatrix>(in, bindMatrix);
      MObject oMatrixData = fnMatrixData.create(bindMatrix, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      plugSampleBindMatrix.elementByLogicalIndex(logicalIndex, &status).setMObject(oMatrixData);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      // Triangle vertices
      MIntArray triangleVertices;
      ReadAttribute<int, MIntArray>(in, triangleVertices);
      MObject oNumericData = fnNumericData.create(MFnNumericData::k3Int, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      status = fnNumericData.setData3Int(triangleVertices[0], triangleVertices[1],
                                         triangleVertices[2]);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      plugTriangleVerts.elementByLogicalIndex(logicalIndex, &status).setMObject(oNumericData);
      CHECK_MSTATUS_AND_RETURN_IT(status);

      // Barycentric coordinates
      MFloatArray coords;
      ReadAttribute<float, MFloatArray>(in, coords);
      oNumericData = fnNumericData.create(MFnNumericData::k3Float, &status);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      status = fnNumericData.setData3Float(coords[0], coords[1], coords[2]);
      CHECK_MSTATUS_AND_RETURN_IT(status);
      plugBarycentricWeights.elementByLogicalIndex(logicalIndex, &status).setMObject(oNumericData);
      CHECK_MSTATUS_AND_RETURN_IT(status);
    }
  }
  MGlobal::displayInfo("Wrap binding imported.");

  return MS::kSuccess;
}
