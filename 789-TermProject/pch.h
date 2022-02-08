#pragma once

























//////////////////////////////////......
#define _CRT_SECURE_NO_DEPRECATE
#include <inttypes.h>
#define HAVE_INT8_T
#include <cmath>
#include <xlocnum>

#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>
#include <cmath>
#include <string>
#include <sstream>
#include <limits>
#include <chrono>
#include <limits>
#include <algorithm>

#include <Eigen3.4-rc1/Eigen>
#include <Eigen3.4-rc1/Geometry>

//#include <Eigen3.2.6/Eigen>
//#include <Eigen3.2.6/Geometry>
//
//#include <Eigen3.0.0/Eigen>
//#include <Eigen3.0.0/Geometry>

typedef unsigned int uint;
typedef long unsigned int luint;



#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoIndexedLineSet.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoTransform.h>

#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>
#include <Inventor/Win/devices/SoWinMouse.h>
#include <Inventor/Win/devices/SoWinKeyboard.h>
#include <Inventor/Win/SoWinCursor.h>

#include <Inventor/Win/SoWinRenderArea.h>
#include <Inventor/Win/devices/SoWinMouse.h>
#include <Inventor/events/SoMouseButtonEvent.h>

using std::string;
using std::vector;
using namespace Eigen;

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

//////////////////////////////////////////
typedef Eigen::Vector3f Vec3; typedef Vec3 vec3;
typedef Eigen::Matrix3f Mat3; typedef Mat3 mat3;
typedef Eigen::Vector4f Vec4; typedef Vec4 vec4;
typedef Eigen::Matrix4f Mat4; typedef Mat4 mat4;
typedef Eigen::Quaternion<float> Quat; typedef Quat quat;

class Mesh; class ARAP; class Painter;
class Vertex; class Edge;class Triangle; 
class Tetrahedron;
////////////////////////////////////////

