#pragma once
#include "pch.h"
//#include "Data.h"


// Decleration

// Definition

class  Vertex
{
public:
	Vector3f coord; // 3d coordinates etc
	Vector3f normal;
	Vector3f color;

	int  idx;   // who am i; verts[idx]

	vector< int > vertList; //adj vertices;
	vector< int > triList;
	vector< int > edgeList;

	// cons
	Vertex(int i,  float x,  float y,  float z);
	Vertex(int i, Vertex* vert);
};

class Edge
{
public:
	int idx; //edges[idx]
	int v1i, v2i; //endpnts
	int t1i, t2i; // right and right tris

	//	float length;   // USING EUCLIDIAN DISTANCE CALCULATION WHEN NECESSARY

	Edge(int id, int v1, int v2);
	Edge(int id, int v1, int v2, int t1, int t2);

	int isNotSingle();
	float getLength(Mesh* msh);
};

class Triangle {
public:
	int idx; //tris[idx]
	int v1i, v2i, v3i;

	Triangle(int id, int v1, int v2, int v3);
	float getArea(Mesh* msh);
};

class Mesh
{
public:
	vector< Vertex* > verts;
	vector< Triangle* > tris;
	vector< Edge* > edges;

	vector< int > samples;
	vector< Edge* > sedges;

	float	minx, maxx, 
			miny, maxy,
			minz, maxz;


	void addTriangle(int v1, int v2, int v3);
	void addEdge(int v1, int v2);
	void addVertex(float x, float y, float z);
	void addVertex(Vec3 nP);
	bool makeVertsNeighbor(int v1i, int v2i);
	void edgeUpdate(int v1, int v2, int id);

	float getLength(int v1, int v2);
	float getArea(Triangle* tri);
	void  loadOff(string name);
	void  loadObj(const char* name);
	void  changeOff(string name);
	void  changeObj(const char* name);
	void CoMToOrigin();
	void  push_left(float x);
	void push_right(float x);

	void emptyMesh();

	Mesh();
	Mesh(string name);
   ~Mesh();
};
