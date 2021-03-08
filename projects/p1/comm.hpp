//===========================================================================//
//                                                                           //
// Copyright(c) 2018 Qi Wu (Wilson)                                          //
// University of California, Davis                                           //
// MIT Licensed                                                              //
//                                                                           //
//===========================================================================//

#include <iostream>
#include <string>
#include <vector>

struct Point {
  float x, y, z;
};

struct Polygon {
  std::vector<Point> points;
  float scale = 1.f;
  float tx = 1.f, ty = 1.f;
  float rotate;
  float cx, cy;
  float TrFloatx = 1.f;//default positions for polygons that can be altered with GUI
  float TrFloaty = 1.f;
  float ScFloat = 1.f;
  float AnFloat = 0.f;
  int size() { return points.size(); }
  float translate_x() { return tx + cx; }
  float translate_y() { return ty + cy; }
  float obj_x(int i) { return points[i].x - cx; }
  float obj_y(int i) { return points[i].y - cy; }
};


extern std::vector<Polygon> polygons;

void ReadFile(const std::string& input);

void ErrorCallback(int error, const char* description);

void MakePix(void* window, int x, int y);

