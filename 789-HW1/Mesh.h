#pragma once

#define _CRT_SECURE_NO_DEPRECATE
#include <inttypes.h>
#define HAVE_INT8_T
#include <cmath>
#include <xlocnum>
#include <Inventor/SbVec3f.h>
#include <iostream>
#include <vector>

using namespace std;

struct Vec3 {
	float x, y, z;
	Vec3() : x(0), y(0), z(0) {};
	Vec3(float x1, float y1, float z1) : x(x1), y(y1), z(z1) {};


  // For sbvec3 addition
	// SbVec3f &operator+(const SbVec3f& rhs) {                 
	// 	return SbVec3f(rhs[0] + x, rhs[1] + y, rhs[2] + z);
	// }

	float * toFloat3(){
		float* array = new float[3];
		
		array[0] = x;
		array[1] = y;
		array[2] = z;
		
		return array;
	}
  
  // Indexing
	float operator[](int idx){
		switch(idx){
		case 0 :
			return x;
		case 1:
			return y;
		case 2:
			return z;
		default:
			return -1;
		}
	}
};

struct Vertex
{
	// float* coords, * normals; //3d coordinates etc
	Vec3 coords,  normals; //3d coordinates etc
	int idx; //who am i; verts[idx]

	vector< int > vertList; //adj vertices;
	vector< int > triList; 
	vector< int > edgeList; 

	SbVec3f color;
	
//	Vertex(int i, float* c) : idx(i), coords(Vec3(c[0] ,c[1], c[2])) {};
	Vertex(int i, float x, float y, float z) : idx(i), coords(Vec3(x,y,z)) {};

};


struct Edge
{
	int idx; //edges[idx]
	int v1i, v2i; //endpnts
//	float length;   // USING EUCLIDIAN DISTANCE CALCULATION WHEN NECESSARY
  
	Edge(int id, int v1, int v2) : idx(id), v1i(v1), v2i(v2) {};
};

struct Triangle
{
	int idx; //tris[idx]
	int v1i, v2i, v3i;
  
	Triangle(int id, int v1, int v2, int v3) : idx(id), v1i(v1), v2i(v2), v3i(v3) {};
};

class Mesh
{
private:
	void addTriangle(int v1, int v2, int v3);
	void addEdge(int v1, int v2);
	void addVertex(float x, float y, float z);
	bool makeVertsNeighbor(int v1i, int v2i);
public:
	vector< Vertex* > verts;
	vector< Triangle* > tris;
	vector< Edge* > edges;

  // painting
	vector< Edge* > sedges;
	vector< int > samples;
	
	Mesh() {} ;
  
  // Utility
	float getLength(int v1,int v2);
	void createCube(float side);
	void loadOff(char* name);
};
