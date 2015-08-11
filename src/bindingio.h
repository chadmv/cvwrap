#ifndef CVWRAP_BindingIO_H
#define CVWRAP_BindingIO_H

#include <maya/MMatrix.h>
#include <maya/MObject.h>
#include <maya/MString.h>
#include <fstream>

/**
  The BindingIO is used to import and export binding information from a wrap node.
*/
class BindingIO {
 public:
  /**
    Exports the binding information to disk.
  */
  MStatus ExportBinding(std::ofstream& out, MObject& oWrapNode);

  /**
    Imports the binding information from disk.
  */
  MStatus ImportBinding(std::ifstream& in, MObject& oWrapNode);

  const static float kWrapFileVersion;
};

/**
  Convenience function to write a Maya array type to a binary stream.
  @param[in] out Output stream.
  @param[in] attribute A Maya array type.
*/
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

/**
  Template specialization for MMatrix because they need to be read differently from
  normal array types.
*/
template <>
void WriteAttribute<double, MMatrix>(std::ofstream &out, const MMatrix& attribute);

/**
  Convenience function to read a Maya array type from a binary stream.
  @param[in] in Input stream.
  @param[out] attribute A Maya array type.
*/
template <typename dataType, typename container>
void ReadAttribute(std::ifstream &in, container &attribute) {
  attribute.clear();
  unsigned int length;
  in.read((char *)(&length), sizeof(length));
  if (length > 0) {
    attribute.setLength(length);
    dataType* pValues = new dataType[length];
    in.read((char *)pValues, length * sizeof(dataType));
    for (unsigned int i = 0; i < length; i++) {
      attribute[i] = pValues[i];
    }
    delete [] pValues;
  }
}

template <>
void ReadAttribute<double, MMatrix>(std::ifstream &in, MMatrix &matrix);


#endif