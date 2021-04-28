#pragma once

#include <iostream>
#include <chrono>
#include <fstream>

#include "Painter.h"
#include "Mesh.h"
#include "Heap.h"
#include "FibHeap.h"

#define INF 2e8
#define UNDF -1

// Q1
#define PATH_START (n-2)
#define PATH_END (n/2)


// Q3
#define SAMP_SIZE 10
#define EXTR_SIZE 15
#define QUALITY 100
#define TOLERANS (0.75)
#define RADIUS 3.0f


class HW1 {
private:
	Mesh* cm;
	int n;
public:
	HW1(Mesh* obje) : cm(obje), n(obje->verts.size()){hw1Main();};

  // Main funtion For HW1
	void hw1Main();

    // VISUALISATION
	void Visualize(); // visualize driver
	void dijVisualize();
	void fpsVisualize();
	void symVisualize();


    //// Q1
    // DIJKSTRA 
	float* aDijkstra(int startingId);     // array implementation
	distpair * hDijkstra(int startingId); // heap  implementation
	float* fhDijkstra(int startingId);    // fiboHeap implementat
	int* fhPathDijkstra(int startingIdx, int finishingIdx); // fibo Path version
	
	// write on file  &  print the time spent for dijkstras (q1)
	void write();
	void timings();

    //// Q2
	int * FPS(); 	// Farthest Point Sampling 
	int AGD();    // Average Geodesiz Distance
	int globalMaxima(float * array);

    /// Q3
	int * sym_inv();
	bool localMaxima(float * array, int id);
	float tolerans(float * array, int id);
};
