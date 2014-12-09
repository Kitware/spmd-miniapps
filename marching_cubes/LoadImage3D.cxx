#include "Image3D.h"
#include "LoadImage3D.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <fstream>
#include <iostream>
#include <sstream>

void loadImage3D(const char *vtiFileName, Image3D *image)
{
  std::ifstream stream;
  stream.open(vtiFileName);

  boost::property_tree::ptree xmlTree;
  boost::property_tree::read_xml(stream, xmlTree);
  stream.close();

  std::stringstream strParseStream;
  std::string value;

  int extent[6], dim[3], npoints;
  double spacing[3], origin[3];

  value = xmlTree.get<std::string>("VTKFile.ImageData.<xmlattr>.WholeExtent");
  strParseStream.str(value);
  strParseStream.seekg(0);
  strParseStream >> extent[0] >> extent[1] >> extent[2] >> extent[3]
                 >> extent[4] >> extent[5];

  dim[0] = extent[1] - extent[0] + 1;
  dim[1] = extent[3] - extent[2] + 1;
  dim[2] = extent[5] - extent[4] + 1;
  npoints = dim[0] * dim[1] * dim[2];

  value = xmlTree.get<std::string>("VTKFile.ImageData.<xmlattr>.Origin");
  strParseStream.str(value);
  strParseStream.seekg(0);
  strParseStream >> origin[0] >> origin[1] >> origin[2];

  value = xmlTree.get<std::string>("VTKFile.ImageData.<xmlattr>.Spacing");
  strParseStream.str(value);
  strParseStream.seekg(0);
  strParseStream >> spacing[0] >> spacing[1] >> spacing[2];

  image->setDimension(dim[0], dim[1], dim[2]);
  image->setOrigin(origin[0], origin[1], origin[2]);
  image->setSpacing(spacing[0], spacing[1], spacing[2]);
  image->allocate();

#if 1
  value =
    xmlTree.get<std::string>("VTKFile.ImageData.Piece.PointData.DataArray");
  strParseStream.str(value);
  strParseStream.seekg(0);
  double *buffer = image->getData();
  for (int i = 0; i < npoints; ++i, ++buffer)
    {
    strParseStream >> *buffer;
    }
#else
  std::vector<short> points(npoints);

  points = xmlTree.get<std::vector<short> >("VTKFile.ImageData.Piece.PointData.DataArray");

  double *buffer = image->getData();
  for (int i = 0; i < npoints; ++i, ++buffer)
    {
    *buffer = points[i];
    }
#endif
}

