//===========================================================================//
//                                                                           //
// Copyright(c) 2018 Qi Wu (Wilson)                                          //
// University of California, Davis                                           //
// MIT Licensed                                                              //
//                                                                           //
//===========================================================================//

#include "util.hpp"
#include "comm.hpp"
#include <stdlib.h>
#include <math.h>
using std::round;
// force imgui to use GLAD
//#define IMGUI_IMPL_OPENGL_LOADER_GLAD
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl2.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;  

static int e = 0;
static float xMinFloat = 0.f;  //These 4 are the display window for world to pixel conversion
static float xMaxFloat = 800.f;
static float yMinFloat = 0.f;
static float yMaxFloat = 600.f;
static float xwMin = 0.f; //These 4 are the world window
static float xwMax = 16.f;
static float ywMin = 0.f;
static float ywMax = 12.f;
static int render = 0;
std::vector<Polygon> polygons;
std::vector<Polygon> outPoly; //vector to output to "input1.txt"
int pid = 0; // current polygon
static bool show_gui = false;
bool buf[800][600]; //boolean frame buffer
class wcPt2D {
  public:
    GLfloat x1, y1;
};
static bool
CapturedByGUI()
{
    ImGuiIO& io = ImGui::GetIO();
    return (io.WantCaptureMouse);
}



void
KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
  //---------------------------------------------------
  // TODO finish the implementation starting from here
  //---------------------------------------------------
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
    // close window
    glfwSetWindowShouldClose(window, GLFW_TRUE);
  }
  else if (key == GLFW_KEY_G && action == GLFW_PRESS) {
    show_gui = !show_gui;
  }
}

void
CursorPositionCallback(GLFWwindow* window, double xpos, double ypos)
{
  //---------------------------------------------------
  // TODO finish the implementation starting from here
  //---------------------------------------------------
  if (!CapturedByGUI()) {
    int left_state  = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
    int right_state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);
    // left click
    if (left_state == GLFW_PRESS) {
    }
    else {
    }
    
    // right click
    if (right_state == GLFW_PRESS) {
    }
    else {
    }
  }
}

static void
WindowSizeCallback(GLFWwindow* window, int width, int height)
{
}

void DDA(void* _window, const int x0, const int y0, const int xEnd, const int yEnd)
{
    int dx = xEnd - x0; //This code is DDA algorithm from the textbook
    int dy = yEnd - y0;
    float xIncrement, yIncrement;
    float x = x0, y = y0;
    int steps;

    if(fabs(dx) > fabs(dy))
    {
        steps = fabs(dx);
    }
    else
    {
        steps = fabs(dy);
    }

    xIncrement = float(dx) / float(steps);
    yIncrement = float(dy) / float(steps);
    MakePix(_window, x, y);
    buf[x0][y0] = true;
    for(int k = 0; k < steps; k++)
    {
        x += xIncrement;
        y += yIncrement;
        MakePix(_window, round(x), round(y));
        if(dy != 0){
        int bufferX = int(round(x));
        int bufferY = int(round(y));
        buf[bufferX][bufferY] = true;
        }
    }
}

void BLD(void* _window, int x0, int y0, int xEnd, int yEnd){
    int dx = fabs(xEnd - x0); //This is the Bresenham algorithm from the book that is modified for special cases.
    int dy = fabs(yEnd - y0);
    int p = 2 * dy - dx;
    int twoDy = 2 * dy, twoDyMinusDx = 2 * (dy - dx);
    int p2 = 2 * dx - dy;
    int twoDx = 2 * dx, twoDxMinusDy = 2 * (dx - dy);
    int x, y;
    int steps;
    x = x0;
    y = y0;
    MakePix(_window, x, y); //every time pixel is made, it gets added to the frame buffer except for horizontal line
    if(dy != 0){
      int bufferX = int(x);
      int bufferY = int(y);
      buf[bufferX][bufferY] = true;
    }
    if(dx > dy){
        steps = dx;
        for(int i = 0; i < steps; i++){
          if(x0 > xEnd){ //for cases where x0 is to the right of xEnd
            x--;
          }
           else{
             x++;
           }
            if(p < 0){ //for negative cases
                p += twoDy;
            }
            else{
              if(y0 > yEnd){
                y--;
              }
              else{
                y++;
              }
                p += twoDyMinusDx;
            }
            MakePix(_window, x, y);
            if(dy != 0){
              int bufferX = int(x);
              int bufferY = int(y);
              buf[bufferX][bufferY] = true;
            }
        }
    }
    else{ //copy of the code above but for slope > 1 cases.
        steps = dy;
        for(int i = 0; i < steps; i++){
          if(y0 > yEnd){
            y--;
          }
          else{
            y++;
          }
          if(p2 < 0){
            p2 += twoDx;
          }
          else{
            if(x0 > xEnd){
              x--;
            }
            else{
              x++;
            }
            p2 += twoDxMinusDy;
            }
            MakePix(_window, x, y);
            if(dy != 0){
            int bufferX = int(x);
            int bufferY = int(y);
            buf[bufferX][bufferY] = true;
            }
        }

    }
}

/* ... */
Point Transformation(double _matrix[3][3], const Point& p) //tranform class that all types of transform will call.
{//takes in 3x3 matrix and multiplied by point p.
    double sum, aPoint[3], newPoint[3];
    aPoint[0] = p.x, aPoint[1] = p.y, aPoint[2] = 1.0;
    for ( int i = 0 ; i < 3 ; i++ ) {
        sum = 0.0 ;
        for ( int j = 0 ; j < 3 ; j++ ) {
            sum += aPoint[j] * _matrix[i][j] ;
        }
        newPoint[i] = sum ;
    }
    Point endPoint;
    endPoint.x = newPoint[0], endPoint.y = newPoint[1], endPoint.z = newPoint[2];
    return endPoint;
}

Point Translate(const Point& p, const Point& t){ //uses point and translate vector to translate
  double matrix[3][3] = {{1.0, 0.0, t.x},{0.0, 1.0, t.y},{0.0,0.0,1.0}} ;
  Point prime = Transformation(matrix, p);
  return prime;
}

Point Scale(const Point& p, float s, float cx, float cy){//scale and rotate by centroid
  double matrix[3][3] = {{s, 0.0, 0.0},{0.0, s, 0.0},{0.0,0.0,1.0}} ;
  Point translate;
  translate.x = cx * -1;
  translate.y = cy * -1;
  Point initPoint = Translate(p, translate);//translate point p to origin
  Point dash = Transformation(matrix, initPoint);//scale it in the origin
  translate.x = cx;
  translate.y = cy;
  Point prime = Translate(dash, translate);//translate back
  return prime;
}

Point Rotate(const Point& p, const float& angle, float cx, float cy){
  double matrix[3][3] = {{cos(angle), sin(angle) * -1, 0.0},{sin(angle), cos(angle), 0.0},{0.0,0.0,1.0}};
  Point translate;
  translate.x = cx * -1;
  translate.y = cy * -1;
  Point initPoint = Translate(p, translate);
  Point dash = Transformation(matrix, initPoint);
  translate.x = cx;
  translate.y = cy;
  Point prime = Translate(dash, translate);
  return prime;
}

Point centerPoint(std::vector<Point> point){//recalculates the center point for every translation
    float cx = 0.f, cy = 0.f;
    Point p;
    for (auto& pt : point) {
      cx += pt.x;
      cy += pt.y;
    }
    cx /= (float)point.size();
    cy /= (float)point.size();
    p.x = cx;
    p.y = cy;
    return p;
}

Point ViewPort(const Point& p){//transfers point p from world to pixel coordinates
  Point v;
  float sx = (800.0)/(xwMax - xwMin);
  float sy = (600.0)/(ywMax - ywMin);
  v.x = 0.0 + (float)((p.x - xwMin) * sx);
  v.y = 0.0 + (float)((p.y - ywMin) * sy);
  return v;
}

const GLint winLeftBitCode = 0x1;
const GLint winRightBitCode = 0x2;
const GLint winBottomBitCode = 0x4;
const GLint winTopBitCode = 0x8;
inline GLint inside (GLint code) { return GLint (!code); }
inline GLint reject (GLint code1, GLint code2)
  { return GLint (code1 & code2); }
inline GLint accept (GLint code1, GLint code2)
  { return GLint (!(code1 | code2)); }
GLubyte encode (wcPt2D pt, wcPt2D winMin, wcPt2D winMax)
{
  GLubyte code = 0x00;
  if (pt.x1 < winMin.x1)
    code = code | winLeftBitCode;
  if (pt.x1 > winMax.x1)
    code = code | winRightBitCode;
  if (pt.y1 < winMin.y1)
    code = code | winBottomBitCode;
  if (pt.y1 > winMax.y1)
    code = code | winTopBitCode;
  return (code);
}
void swapPts (wcPt2D * p1, wcPt2D * p2)
{
  wcPt2D tmp;
  tmp = *p1; *p1 = *p2; *p2 = tmp;
}
void swapCodes (GLubyte * c1, GLubyte * c2)
{
  GLubyte tmp;
  tmp = *c1; *c1 = *c2; *c2 = tmp;
}
std::vector<Point> lineClipCohSuth (wcPt2D winMin, wcPt2D winMax, wcPt2D p1, wcPt2D p2)
{//Cohen-Sutherland line clipping from textbook but modified to return a vector of points for line drawing
  GLubyte code1, code2;
  GLint done = false, plotLine = false;
  GLfloat m;
  while (!done) {
    code1 = encode (p1, winMin, winMax);
    code2 = encode (p2, winMin, winMax);
    if (accept (code1, code2)) {
      done = true;
      plotLine = true;
    }
    else
      if (reject (code1, code2))
        done = true;
      else {
/* Label the endpoint outside the display window as p1. */
        if (inside (code1)) {
          swapPts (&p1, &p2);
          swapCodes (&code1, &code2);
        }
/* Use slope m to find line-clipEdge intersection. */
        if (p2.x1 != p1.x1)
          m = (p2.y1 - p1.y1) / (p2.x1 - p1.x1);
        if (code1 & winLeftBitCode) {
          p1.y1 += (winMin.x1 - p1.x1) * m;
          p1.x1 = winMin.x1;
        }
        else
          if (code1 & winRightBitCode) {
            p1.y1 += (winMax.x1 - p1.x1) * m;
            p1.x1 = winMax.x1;
        }
        else
          if (code1 & winBottomBitCode) {
/* Need to update p1.x for nonvertical lines only. */
            if (p2.x1 != p1.x1)
              p1.x1 += (winMin.y1 - p1.y1) / m;
              p1.y1 = winMin.y1;
          }
          else
            if (code1 & winTopBitCode) {
              if (p2.x1 != p1.x1)
                p1.x1 += (winMax.y1 - p1.y1) / m;
                p1.y1 = winMax.y1;
          }

}
}
  std::vector<Point> points;
  Point pOne;
  Point pTwo;
  if(!reject(code1, code2)){
    pOne.x = p1.x1;
    pOne.y = p1.y1;
    pTwo.x = p2.x1;
    pTwo.y = p2.y1;
    points.push_back(pOne);
    points.push_back(pTwo);
  }
  return points;
}

typedef enum { Left, Right, Bottom, Top } Boundary;
const GLint nClip = 4;
const Boundary Plus [] =
{
  Left, Right, Bottom, Top
};
GLint inside (wcPt2D p, Boundary b, wcPt2D wMin, wcPt2D wMax)
{
  switch (b) {
    case Left: if (p.x1 < wMin.x1) return (false); break;
    case Right: if (p.x1 > wMax.x1) return (false); break;
    case Bottom: if (p.y1 < wMin.y1) return (false); break;
    case Top: if (p.y1 > wMax.y1) return (false); break;
}
  return (true);
}
GLint cross (wcPt2D p1, wcPt2D p2, Boundary winEdge, wcPt2D wMin, wcPt2D wMax)
{
  if (inside (p1, winEdge, wMin, wMax) == inside (p2, winEdge, wMin, wMax))
    return (false);
  else return (true);
}
wcPt2D intersect (wcPt2D p1, wcPt2D p2, Boundary winEdge,
wcPt2D wMin, wcPt2D wMax)
{
  wcPt2D iPt;
  GLfloat m;
  if (p1.x1 != p2.x1) m = (p1.y1 - p2.y1) / (p1.x1 - p2.x1);
  switch (winEdge) {
  case Left:
  iPt.x1 = wMin.x1;
  iPt.y1 = p2.y1 + (wMin.x1 - p2.x1) * m;
  break;
  case Right:
  iPt.x1 = wMax.x1;
  iPt.y1 = p2.y1 + (wMax.x1 - p2.x1) * m;
  break;
  case Bottom:
  iPt.y1 = wMin.y1;
  if (p1.x1 != p2.x1) iPt.x1 = p2.x1 + (wMin.y1 - p2.y1) / m;
  else iPt.x1 = p2.x1;
  break;
  case Top:
  iPt.y1 = wMax.y1;
  if (p1.x1 != p2.x1) iPt.x1 = p2.x1 + (wMax.y1 - p2.y1) / m;
  else iPt.x1 = p2.x1;
  break;
}
return iPt;
}
void clipPoint (wcPt2D p, Boundary winEdge, wcPt2D wMin, wcPt2D wMax,
wcPt2D * pOut, int * cnt, wcPt2D * first[], wcPt2D * s)
{
  wcPt2D iPt;
/* If no previous point exists for this clipping boundary,
* save this point.
*/
  if (!first[winEdge])
    first[winEdge] = &p;
  else
/* Previous point exists. If p and previous point cross
* this clipping boundary, find intersection. Clip against
* next boundary, if any. If no more clip boundaries, add
* intersection to output list.
*/
    if (cross (p, s[winEdge], winEdge, wMin, wMax)) {
      iPt = intersect (p, s[winEdge], winEdge, wMin, wMax);
      if (winEdge < Top)
          clipPoint (iPt, Plus[winEdge + 1], wMin, wMax, pOut, cnt, first, s);
      else {
        pOut[*cnt] = iPt; (*cnt)++;
}
}
/* Save p as most recent point for this clip boundary. */
    s[winEdge] = p;
/* For all, if point inside, proceed to next boundary, if any. */
    if (inside (p, winEdge, wMin, wMax))
      if (winEdge < Top)
        clipPoint (p, Plus[winEdge + 1], wMin, wMax, pOut, cnt, first, s);
      else {
        pOut[*cnt] = p; (*cnt)++;
}
}
void closeClip (wcPt2D wMin, wcPt2D wMax, wcPt2D * pOut,
GLint * cnt, wcPt2D * first [ ], wcPt2D * s)
{
  wcPt2D pt;
  Boundary winEdge;
  for (int i = 0; i < Top; i++) {
    if (cross (s[winEdge], *first[winEdge], winEdge, wMin, wMax)) {
      pt = intersect (s[winEdge], *first[winEdge], winEdge, wMin, wMax);
      if (winEdge < Top)
        clipPoint (pt, Plus[i+1], wMin, wMax, pOut, cnt, first, s);
      else {
        pOut[*cnt] = pt; (*cnt)++;
}
}
}
}


GLint polygonClipSuthHodg (wcPt2D wMin, wcPt2D wMax, GLint n, wcPt2D * pIn, wcPt2D * pOut)
{//Completed version of Sutherland-Hodgeman from the textbook
wcPt2D * first[nClip] = { 0, 0, 0, 0 }, s[nClip];
GLint k, cnt = 0;
for (k = 0; k < n; k++)
  clipPoint (pIn[k], Left, wMin, wMax, pOut, &cnt, first, s);
  closeClip (wMin, wMax, pOut, &cnt, first, s);
  return (cnt);
}


struct edgeBucket {//edge struct for an edge table.
    int yMax;
    int yMin;
    int xMin;
    int xMax;
};

std::vector<edgeBucket> CreateEdges(std::vector<Point> points){//takes in a polygon to construct an edge table
  std::vector<edgeBucket> edgeTable;
  int yMax; int yMin; int xMin; int xMax;
  edgeBucket edge;
  for(int i = 1; i < points.size(); i++){//for loop generates edges for every point
    int dx = points[i].x - points[i-1].x;
    int dy = points[i].y - points[i-1].y;
    if(points[i].y > points[i-1].y){
      yMax = points[i].y;
      yMin = points[i-1].y;
      xMin = points[i-1].x;
      xMax = points[i].x;
    }
    else{
      yMax = points[i-1].y;
      yMin = points[i].y;
      xMin = points[i].x;
      xMax = points[i-1].x;
    }
    if(dy != 0){ //non horizontal points are included in the edge table
      edge.yMax = yMax;
      edge.yMin = yMin;
      edge.xMin = xMin;
      edge.xMax = xMax;
      edgeTable.push_back(edge);
    }
  }
  int dx = points[points.size()-1].x - points[0].x;//makes edge for first and last points
  int dy = points[points.size()-1].y - points[0].y;
  if(points[points.size()-1].y > points[0].y){
    yMax = points[points.size()-1].y;
    yMin = points[0].y;
    xMin = points[0].x;
    xMax = points[points.size()-1].x;
  }
  else{
    yMax = points[0].y;
    yMin = points[points.size()-1].y;
    xMin = points[points.size()-1].x;
    xMax = points[0].x;
  }
  if(dy != 0){
    edge.yMax = yMax;
    edge.yMin = yMin;
    edge.xMin = xMin;
    edge.xMax = xMax;
    edgeTable.push_back(edge);
  }
  return (edgeTable);
}
void Render(void* _window, Point bMin, Point bMax, std::vector<Point> shape){//Rasterizes based on boolean frame buffer.
  int parity = -1; //takes in bounding box and list of vertices
  int parity2 = -1;//1 means inside//-1 means outside
  for(int j = bMin.y; j <= bMax.y; j++){//Next two for loops are for clipping cases
    if(xMinFloat > bMin.x){
      if(buf[int(xMinFloat)][j] == true){
        parity *= -1;
      }
      if(parity == 1){
        buf[int(xMinFloat)][j] = true;
      }
    }
    if(xMaxFloat < bMax.x){
      if(buf[int(xMaxFloat)][j] == true){
        parity2 *= -1;
      }
      if(parity2 == 1){
        buf[int(xMaxFloat)][j] = true;
      }
    }
  }
  parity = -1;
  parity2 = -1;
  for(int i = bMin.x; i <= bMax.x; i++){
    if(yMinFloat > bMin.y){
      if(buf[i][int(yMinFloat)] == true){
        parity *= -1;
      }
      if(parity == 1){
        buf[i][int(yMinFloat)] = true;
      }
    }
    if(yMaxFloat < bMax.y){
      if(buf[i][int(yMaxFloat)] == true){
        parity2 *= -1;
      }
      if(parity2 == 1){
        buf[i][int(yMaxFloat)] = true;
      }
    }
  }
  std::vector<edgeBucket> edgeTable = CreateEdges(shape);//the purpose of edge table is to detect vertices.
  for(int j = bMin.y; j <= bMax.y; j++){ //j is the scan line
    int parity = -1;
    for(int i = bMin.x; i <= bMax.x; i++){
      if(buf[i][j] == true){ //whenever scanline touches the intersection
        parity *= -1;
        if(buf[i+1][j] == true){ //for horizontal cases
          parity *= -1;
        }
        int verticesMin = 0;
        int verticesMax = 0;
        for(int f = 0; f < edgeTable.size(); f++){ //for vertices cases
          if(edgeTable[f].xMin == i && edgeTable[f].yMin == j){
            verticesMin++;
        }
          if(edgeTable[f].xMax == i && edgeTable[f].yMax == j){
            verticesMax++;
        }
        }
        if(verticesMin == 2){
          parity *= -1;
          }
        else if(verticesMax == 2){
          parity *= -1;
        }
      }
      if(parity == 1){ //makes pixels if inside polygon and inside viewport
        if(i >= xMinFloat && i <= xMaxFloat && j >= yMinFloat && j <= yMaxFloat){
            MakePix(_window, i, j);
        }
      }
        buf[i][j] = false; //resets the buffer
      }
    }
  }

void
DrawCall(GLFWwindow* window)
{
  //---------------------------------------------------
  // TODO finish the implementation starting from here
  //---------------------------------------------------
  Point prev;
  Point big;
  int track = 0;
  prev.x = 0;
  prev.y = 0;
  wcPt2D worldMin; wcPt2D worldMax;
  worldMax.x1 = xMaxFloat;
  worldMax.y1 = yMaxFloat;
  worldMin.x1 = xMinFloat;
  worldMin.y1 = yMinFloat;
  wcPt2D p1; wcPt2D p2;
  for (auto& p : polygons){
    Polygon out;
    Point boxMin; boxMin.x = INFINITY; boxMin.y = INFINITY;
    Point boxMax; boxMax.x = -INFINITY; boxMax.y = -INFINITY;
    std::vector<Point> shape;
    prev.x = p.TrFloatx;
    prev.y = p.TrFloaty;
    for (int i = 0; i < p.points.size(); i++){//scales every point
      big = Scale(p.points[i], p.ScFloat, p.cx, p.cy);
      shape.push_back(big);
    }
    for (int i = 0; i < shape.size(); i++){//rotates every point
      big = Rotate(shape[i], p.AnFloat, p.cx, p.cy);
      shape[i] = big;
    }
    for (int i = 0; i < shape.size(); i++){
      big = Translate(shape[i], prev);//translates every point
      shape[i] = big;
      out.points.push_back(shape[i]);//adds to output vector
      Point view = ViewPort(shape[i]); //converts world points to pixels
      shape[i] = view;
      if(shape[i].x < boxMin.x){//boxMin and boxMax creates bounding box for rasterizing
        boxMin.x = shape[i].x;
      }
      if(shape[i].y < boxMin.y){
        boxMin.y = shape[i].y;
      }
      if(shape[i].x > boxMax.x){
        boxMax.x = shape[i].x;
      }
      if(shape[i].y > boxMax.y){
        boxMax.y = shape[i].y;
      }
      if(i > 0){
        p1.x1 = shape[i-1].x;//clips the lines
        p1.y1 = shape[i-1].y;
        p2.x1 = shape[i].x;
        p2.y1 = shape[i].y;
        std::vector<Point> get = lineClipCohSuth(worldMin, worldMax, p1, p2);
        if(get.size() != 0){
          if(e == 0){//draws the lines
            DDA(window, get[0].x, get[0].y, get[1].x, get[1].y);
          }
          else{
            BLD(window, get[0].x, get[0].y, get[1].x, get[1].y);
          }
        }
      }
    }
    p1.x1 = shape[0].x;//first and last point cases
    p1.y1 = shape[0].y;
    p2.x1 = shape[shape.size()-1].x;
    p2.y1 = shape[shape.size()-1].y;
    std::vector<Point> get = lineClipCohSuth(worldMin, worldMax, p1, p2);
    if(get.size() != 0){
      if(e==0){
        DDA(window, get[0].x, get[0].y, get[1].x, get[1].y);
      }
      else{
        BLD(window, get[0].x, get[0].y, get[1].x, get[1].y);
      }
    }
    if(render == 0){//rasterizes
      Render(window, boxMin, boxMax, shape);
    }
    outPoly[track] = out;//saves all changes to output vector
    track++;
  }
}


int
main(const int argc, const char** argv)
{
  if (argc < 2) 
    throw std::runtime_error("missing input file"); 
  
  ReadFile(argv[1]);
  // Compute the Center of Mass
  for (auto& p : polygons) {        
    float cx = 0.f, cy = 0.f;
    for (auto& pt : p.points) {
      cx += pt.x;
      cy += pt.y;
    }
    cx /= (float)p.points.size();
    cy /= (float)p.points.size();
    p.cx = cx;
    p.cy = cy;
  }
  outPoly.resize(polygons.size());
  //---------------------------------------------------
  // TODO finish the implementation starting from here
  //---------------------------------------------------
  int width = 800, height = 600;
  for(int i = 0; i < 800; i++){ //initializes frame buffer
    for(int j = 0; j < 600; j++){
      buf[i][j] = false;
    }
  }
  // Initialize GLFW
  glfwSetErrorCallback(ErrorCallback);
  if (!glfwInit()) {
    exit(EXIT_FAILURE);
  }

  // Provide Window Hint
  glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

  // OpenGL Setup
  GLFWwindow* window = NULL;
  window = glfwCreateWindow(width, height,
			    "ECS 175 Renderer",
			    NULL, NULL);
  if (!window) {
    glfwTerminate();
    throw std::runtime_error("Failed to create GLFW window");
  }

  // Callback
  glfwSetKeyCallback(window, KeyCallback);
  glfwSetWindowSizeCallback(window, WindowSizeCallback);
  glfwSetCursorPosCallback(window, CursorPositionCallback);

  // Ready
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);

  // Load GLAD symbols
  int err = gladLoadGLLoader((GLADloadproc)glfwGetProcAddress) == 0;
  if (err) {
    throw std::runtime_error("Failed to initialize OpenGL loader!");
  }

  // ImGui
  {
    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    
    // Setup Dear ImGui style
    ImGui::StyleColorsDark(); // or ImGui::StyleColorsClassic();
    
    // Initialize Dear ImGui
    ImGui_ImplGlfw_InitForOpenGL(window, false);
    ImGui_ImplOpenGL2_Init();
  }
  // Execute
  while (!glfwWindowShouldClose(window)) {
    glClear(GL_COLOR_BUFFER_BIT);
    // Draw Background
    DrawCall(window);
    
    // Draw GUI
    {
      // Initialization
      ImGui_ImplOpenGL2_NewFrame();
      ImGui_ImplGlfw_NewFrame();
      ImGui::NewFrame();
      // - Uncomment below to show ImGui demo window
      ImGuiTabBarFlags tab_bar_flags = ImGuiTabBarFlags_None;
      if (ImGui::BeginTabBar("MyTabBar", tab_bar_flags))
      {
        for(int i = 0; i<polygons.size(); i++){//creates GUI tab bar for every polygon
          std::vector<char> arr;
          char a = '0' + i;
          arr.push_back(a);
          const char *b = &a;
        if (ImGui::BeginTabItem(b))
          {
            ImGui::Text("Polygon ID");
            ImGui::Text("%d\n", i);
            ImGui::Text("Line Algorithm");
            ImGui::RadioButton("DDA", &e, 0); ImGui::SameLine();
            ImGui::RadioButton("Bresenham", &e, 1);
            ImGui::Text("Transformation");
            ImGui::SliderFloat("Translation X", &polygons[i].TrFloatx, 0.f, 8.f);
            ImGui::SliderFloat("Translation Y", &polygons[i].TrFloaty, 0.f, 6.f);
            ImGui::SliderFloat("Scaling Factor", &polygons[i].ScFloat, 0.f, 4.f);
            ImGui::SliderFloat("Rotation Angle", &polygons[i].AnFloat, -360.f, 360.f);
            ImGui::Text("Viewport");
            ImGui::SliderFloat("X-Extention Min", &xMinFloat, 0.f, 800.f);
            ImGui::SliderFloat("X-Extension Max", &xMaxFloat, 0.f, 800.f);
            ImGui::SliderFloat("Y-Extention Min", &yMinFloat, 0.f, 600.f);
            ImGui::SliderFloat("Y-Extension Max", &yMaxFloat, 0.f, 600.f);
            ImGui::RadioButton("Render", &render, 0); ImGui::SameLine();
            ImGui::RadioButton("NoRender", &render, 1);
            if(ImGui::Button("Update Input")){
                std::ofstream outFile; //writes to input file when windows close
                outFile.open(argv[1]);
                if(outFile.is_open()){
                  outFile <<outPoly.size();
                  outFile <<"\n";
                  for(auto &p : outPoly){
                    outFile <<"\n";
                    outFile <<p.points.size();
                    outFile <<"\n";
                    for(int i = 0; i < p.points.size(); i++){
                      outFile <<p.points[i].x;
                      outFile <<" ";
                      outFile <<p.points[i].y;
                      outFile <<"\n";
                    }
                  }
                  outFile.close();
                }
                else{
                  std::cout <<"Cannot Open File";
                }
            }
            ImGui::EndTabItem();
          }
        }
            ImGui::EndTabBar();
}
      // Render GUI
      ImGui::Render();
      ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
    }
    glfwSwapBuffers(window);
    glfwPollEvents();
  }

  // Exit
  ImGui_ImplOpenGL2_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();
  
  glfwDestroyWindow(window);
  glfwTerminate();
  return EXIT_SUCCESS;
}
