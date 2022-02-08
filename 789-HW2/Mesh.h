#pragma once

#define _CRT_SECURE_NO_DEPRECATE
#include <inttypes.h>
#define HAVE_INT8_T
#include <cmath>
#include <xlocnum>
#include <Inventor/SbVec3f.h>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;



struct Vec3 {
	float x, y, z;
	Vec3() : x(0), y(0), z(0) {};
	Vec3(float x1, float y1, float z1) : x(x1), y(y1), z(z1) {};

  // For sbvec3 addition
	// SbVec3f &operator+(const SbVec3f& rhs) {                 
	// 	return SbVec3f(rhs[0] + x, rhs[1] + y, rhs[2] + z);
	// }


	
	Vec3 normalize(){
		float mag = sqrt(x*x + y*y + z*z);
		return Vec3(x/mag, y/mag, z/mag);
	}

	float magnitude(){
		return sqrt(x*x + y*y + z*z);
	}
	
	float * toFloat3(){
		float* array = new float[3];
		
		array[0] = x;
		array[1] = y;
		array[2] = z;
		
		return array;
	}
  
  // OPERATORS
	bool operator==(const Vec3 &rhs) { // Scaler Div 
		return (x == rhs.x && y == rhs.y && z == rhs.z);
	}
	Vec3 operator/(float n ){ // Scaler Div 
		return Vec3(x/n,y/n,z/n);
	}
	Vec3 operator^(Vec3 rhs){ // Cross
		return Vec3((y*rhs.z)-(z*rhs.y),  (z*rhs.x)-(x*rhs.z),  (x*rhs.y)-(y*rhs.x));
	}
	Vec3 operator^(float n){ // Scaler Power
		return Vec3(pow(x,n),  pow(y,n), pow(z,n));
	}
	float operator*(Vec3 rhs){ // Dot
		return ((x*rhs.x) + (y*rhs.y) + (z*rhs.z));
	}
	Vec3 operator*(float n ){ // Scalar Prod
		return Vec3(x*n,y*n,z*n);
	}
	Vec3 operator-(const Vec3 &rhs){
		return Vec3(x - rhs.x, y-rhs.y, z-rhs.z);
	}
	Vec3 operator+(const Vec3 &rhs){
		return Vec3(x + rhs.x, y+rhs.y, z+rhs.z);
	}
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
struct Vec3i {
	int x, y, z;
	Vec3i() : x(0), y(0), z(0) {};
	Vec3i(int x1, int y1, int z1) : x(x1), y(y1), z(z1) {};
	Vec3i(float x1, float y1, float z1) : x((int)x1), y((int)y1), z((int)z1) {};

	// For sbvec3 addition
	  // SbVec3if &operator+(const SbVec3if& rhs) {                 
	  // 	return SbVec3if(rhs[0] + x, rhs[1] + y, rhs[2] + z);
	  // }

	Vec3 normalize() {
		float mag = sqrt(x * x + y * y + z * z);
		return Vec3(x / mag, y / mag, z / mag);
	}

	float magnitude() {
		return sqrt(x * x + y * y + z * z);
	}

	float* toFloat3() {
		float* array = new float[3];

		array[0] = x;
		array[1] = y;
		array[2] = z;

		return array;
	}

	// OPERATORS
	Vec3 operator/(float n) { // Scaler Div 
		return Vec3(x / n, y / n, z / n);
	}
	Vec3i operator^(Vec3i rhs) { // Cross
		return Vec3i((y * rhs.z) - (z * rhs.y), (z * rhs.x) - (x * rhs.z), (x * rhs.y) - (y * rhs.x));
	}
	float operator*(Vec3i rhs) { // Dot
		return ((x * rhs.x) + (y * rhs.y) + (z * rhs.z));
	}
	Vec3i operator*(float n) { // Scalar Prod
		return Vec3i(x * n, y * n, z * n);
	}
	Vec3i operator-(Vec3i& rhs) {
		return Vec3i(x - rhs.x, y - rhs.y, z - rhs.z);
	}
	Vec3i operator+(const Vec3i& rhs) {
		return Vec3i(x + rhs.x, y + rhs.y, z + rhs.z);
	}
	float operator[](int idx) {
		switch (idx) {
		case 0:
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

struct Line {
	float m, k;

	Line() : m(0), k(0) {};
	Line(float m1, float k1) : m(m1), k(k1) {};
	Line(Vec3 p1, Vec3 p2) : m((p2.y - p1.y) / (p2.x - p1.x)) { k = p1.y - p1.x * m; };


	bool isOnRight(Vec3& rhs) {
		return (m < 0) ? !(rhs.x * m + k < rhs.y) : (rhs.x * m + k < rhs.y);
	}

};
struct Vertex
{
	// float* coords, * normals; //3d coordinates etc
	Vec3 coords,  normal; //3d coordinates etc
	int idx; //who am i; verts[idx]

	vector< int > vertList; //adj vertices;
	vector< int > triList; 
	vector< int > edgeList; 

	SbVec3f color;
	
//	Vertex(int i, float* c) : idx(i), coords(Vec3(c[0] ,c[1], c[2])) {};
	Vertex(int i, float x, float y, float z) : idx(i), coords(Vec3(x,y,z)), normal(Vec3(0,0,0)), color(SbVec3f(0, 255, 0)) {};
	Vertex(int i, Vertex* vert) : idx(i), coords(Vec3(vert->coords)), normal(Vec3(vert->normal)), color(SbVec3f(0, 255, 0)) {};

};


struct Edge
{
	int idx; //edges[idx]
	int v1i, v2i; //endpnts
	int t1i, t2i; // right and right tris

//	float length;   // USING EUCLIDIAN DISTANCE CALCULATION WHEN NECESSARY

	Edge(int id, int v1, int v2) : idx(id), v1i(v1), v2i(v2), t1i(-1), t2i(-1) {};
	Edge(int id, int v1, int v2, int t1, int t2) : idx(id), v1i(v1), v2i(v2), t1i(t1), t2i(t2) {};

	int isNotSingle(){
		return (t1i == -1)*2 + (t2i == -1)*1;
	}
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
public:
	//private
	void addTriangle(int v1, int v2, int v3);
	void addEdge(int v1, int v2);
	void addVertex(float x, float y, float z);
	void addVertex(Vec3 nP);
	bool makeVertsNeighbor(int v1i, int v2i);
	void edgeUpdate(int v1, int v2, int id);

	//public
	vector< Vertex* > verts;
	vector< Triangle* > tris;
	vector< Edge* > edges;


  // painting
	vector< Edge* > sedges;
	vector< int > samples;
	
	Mesh() {} ;
	~Mesh(){
		for (int i = 0; i < verts.size(); ++i) {
			delete verts[i];
		}
		for (int i = 0; i < edges.size(); ++i) {
			delete edges[i];
		}
		for (int i = 0; i < tris.size(); ++i) {
			delete tris[i];
		}
		for (int i = 0; i < sedges.size(); ++i) {
			delete sedges[i];
		}
	};
  // Utility
	float getArea(Triangle* tri);
	float getLength(int v1,int v2);
	void createCube(float side);
	void loadOff(char* name);
	Vec3 coord(int ix);	

};
