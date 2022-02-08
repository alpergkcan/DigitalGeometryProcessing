#pragma once
#include "Settings.h"
#if HW_2 || FINAL

#include <iostream>
#include <chrono>
#include <fstream>

#include "Painter.h"
#include "Mesh.h"
#include "myHeap.h"
#include "myFibHeap.h"

class HW2 {
public:

// FIELDS
    // private
	bool exit;
	bool paint;
	float alpha;
	bool astar;
	bool subdiv;
	Mesh* cm;
	Mesh* nm;	
	int n;
	
    // public
	HW2(Mesh* obje) : cm(obje), n(obje->verts.size()), exit(false), alpha(ALPHA), paint(false),astar(false), subdiv(false) { hw2Main(); };
//..



// FUNCTIONS
	
	// Main funtion For HW1
	void hw2Main();
	// Main Utilities
	void compare();
	double pathLength(int* path, int s, int t );	
	double surfArea(Mesh* mesh);
	int* fhPathDijkstra(int startingIdx, int finishingIdx);
	Vec3 chooseBary(Triangle* tri, float div = 0.1);
	//// Q1
	void aStar(int s, int t);
	void aStarCross(int s, int t);
	float h(int s, int t);

	//// Q2
	// Subdivisions
	Mesh * _3Subdiv();
	Mesh * _4Subdiv();
	Mesh * _PSubdiv();	
	Mesh * __PSubdiv();

	// utility
	Vec3 findNormal(int v1, int v2, int v3);
	vector<int> getEdges(Triangle* tri);
	vector<int> getEdgesM(Mesh* mesh, Triangle* tri);

	Vec3 projTangent(Vec3 q, Vertex * v);
	Vec3 p(Vec3 b, Vertex** v);
	Vec3 p_star(Vec3 b, Triangle* t);
	
	void triangulate(Mesh* sub, vector<Vertex*> pset);
	vector<Vertex*> triInterior( vector<Vertex *> pset, Vec3 v1, Vec3 v2, Vec3 v3 );
	vector<int> recurTri(Mesh* mesh, vector<Vertex*> pset, Vertex* v1, Vertex* v2, Vertex* v3, vector<bool>& ishandled);
	

	void flip(Mesh* mesh,
		  Triangle* t1, Triangle* t2,
		  Edge * e,
		  int v1, int v2, int v3, int v4);
	
	float minAngleForFlip(Vec3 v1, Vec3 v2, Vec3 v3, Vec3 v4);
	void delanuay(Mesh* mesh, vector<Edge*> edg);

	vector<Vec3i> polygonTriangulation(Mesh* mesh, vector<int> convexHull);
	vector<Vec3> createPointSet(Triangle* tri, int beta, vector<int>pMods); // vector<Vertex*>
	vector<int> grahamScan(vector<Vertex *> pset);
	int parallel(vector<Vertex*> pset);

	Vec3 tangentPlaneCorrespondence(Vec3 p, Vec3 n, Vec3 c);
	float angleInDeg(Vec3 v1, Vec3 v2);
	double angleInGrad(Vec3 v1, Vec3 v2);
	void copyVertex(Mesh *msh, Mesh* sub, int idx);
	bool isInside(Vec3 point, Vec3 v1, Vec3 v2, Vec3 v3);
	vector<int> uniqeSet(vector<int> set );


	vector<Vec3>  basicPointSet(int mod);
	void basicTriangulate(Mesh* sub, vector<Vertex*> pset, int mod);
};

//,,
#endif
