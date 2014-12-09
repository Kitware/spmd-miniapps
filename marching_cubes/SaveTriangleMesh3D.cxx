#include "SaveTriangleMesh3D.h"
#include "TriangleMesh3D.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <fstream>
#include <iostream>
#include <sstream>

void saveTriangleMesh3D(const TriangleMesh3D &mesh, const char *vtpFileName)
{
  boost::property_tree::ptree vtkfile;
  vtkfile.put("VTKFile.<xmlattr>.type", "PolyData");
  vtkfile.put("VTKFile.<xmlattr>.version", "0.1");
  vtkfile.put("VTKFile.<xmlattr>.byte_order", "LittleEndian");

  int nverts = mesh.vertices.size();
  int ninds = mesh.indexes.size();
  int ntris = ninds/3;

  boost::property_tree::ptree polydata;
  polydata.put("Piece.<xmlattr>.NumberOfPoints", nverts);
  polydata.put("Piece.<xmlattr>.NumberOfVerts", 0);
  polydata.put("Piece.<xmlattr>.NumberOfLines", 0);
  polydata.put("Piece.<xmlattr>.NumberOfStrips", 0);
  polydata.put("Piece.<xmlattr>.NumberOfPolys", ntris);

  std::stringstream formatStream;

  boost::property_tree::ptree points;
  points.put("DataArray.<xmlattr>.type", "Float64");
  points.put("DataArray.<xmlattr>.NumberOfComponents", 3);
  points.put("DataArray.<xmlattr>.format", "ascii");
  for (int i = 0; i < nverts; ++i)
    {
    formatStream << mesh.vertices[i].pos[0] << " "
                 << mesh.vertices[i].pos[1] << " "
                 << mesh.vertices[i].pos[2] << "\n";
    }
  points.put("DataArray", formatStream.str());
  formatStream.str(std::string());

  boost::property_tree::ptree pointdata;
  pointdata.put("<xmlattr>.Normals", "normals");
  pointdata.put("DataArray.<xmlattr>.type", "Float64");
  pointdata.put("DataArray.<xmlattr>.NumberOfComponents", 3);
  pointdata.put("DataArray.<xmlattr>.Name", "normals");
  pointdata.put("DataArray.<xmlattr>.format", "ascii");
  for (int i = 0; i < nverts; ++i)
    {
    formatStream << mesh.vertices[i].norm[0] << " "
                 << mesh.vertices[i].norm[1] << " "
                 << mesh.vertices[i].norm[2] << "\n";
    }
  pointdata.put("DataArray", formatStream.str());
  formatStream.str(std::string());

  boost::property_tree::ptree polys;
  boost::property_tree::ptree connectivity, offsets;
  connectivity.put("<xmlattr>.type", "Int32");
  connectivity.put("<xmlattr>.Name", "connectivity");
  connectivity.put("<xmlattr>.format", "ascii");
  for (int i = 0; i < ninds; ++i)
    {
    formatStream << mesh.indexes[i] << "\n";
    }
  connectivity.put_value(formatStream.str());
  formatStream.str(std::string());
  offsets.put("<xmlattr>.type", "Int32");
  offsets.put("<xmlattr>.Name", "offsets");
  offsets.put("<xmlattr>.format", "ascii");
  for (int i = 0, o = 3; i <= ntris; ++i, o += 3)
    {
    formatStream << o << "\n";
    }
  offsets.put_value(formatStream.str());
  formatStream.str(std::string());
  polys.add_child("DataArray", connectivity);
  polys.add_child("DataArray", offsets);

  polydata.add_child("Piece.Points", points);
  polydata.add_child("Piece.PointData", pointdata);
  polydata.add_child("Piece.Polys", polys);
  vtkfile.add_child("VTKFile.PolyData", polydata);

  std::ofstream stream;
  stream.open(vtpFileName);
  boost::property_tree::write_xml(stream, vtkfile,
    boost::property_tree::xml_writer_make_settings('\t', 1));
  stream.close();
}

