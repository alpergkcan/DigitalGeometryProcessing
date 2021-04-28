#include "HW1.h"


void HW1::hw1Main(){
	std::cout << "HOMEWORK ASSIGNMENT 1 " << std::endl << std::endl;
	
choice:	std::cout << "1)  Geodesic on Meshes - Dijkstra  " << std::endl << std::endl;
	std::cout << "2)  Sampling on Meshes - Average Geodesic Distances " << std::endl << std::endl;
	std::cout << "3)  Symmetry-Invariant Sampling " << std::endl << std::endl;
	std::cout << "4)  All At Once.. " << std::endl << std::endl;

	int a;
	std::cin >> a;
	switch(a){
	case 1:
		dijVisualize();
		std::cout << std::endl;  
		write();
		std::cout << std::endl;
		timings();
		std::cout << std::endl<<std::endl;
		break;    
	case 2:
		fpsVisualize();
		std::cout << std::endl<<std::endl;
		break;
	case 3:
		symVisualize();
		break;
	case 4:
		std::cout << std::endl<< std::endl;
		std::cout << "1st Part)  Geodesic on Meshes - Dijkstra  " << std::endl << std::endl;
		dijVisualize();
		std::cout << std::endl;  
		write();
		std::cout << std::endl;
		timings();
		std::cout << std::endl<<std::endl;
		std::cout << "2nd Part)  Sampling on Meshes - Average Geodesic Distances " << std::endl << std::endl;
		fpsVisualize();
		std::cout << std::endl<<std::endl;
		std::cout << "3rd Part)  Symmetry-Invariant Sampling " << std::endl << std::endl;
		symVisualize();
		break;
	default:
		std::cout << "Please Choose Another option" << std::endl;
		goto choice;
	};


}



int HW1::AGD(){
	float * agd = new float[n];
	float * mgd = new float[n];

	for (int i = 0; i < n; ++i){
		agd[i] = 0;
		mgd[i] = INF;
	}
	
	// make QUALITY # of samples  to calculate agd for each v
	int seed = 0;

	float * shortestDistance;
	std::cout << "AGD samples: " << std::endl;	
	for (int i = 0; i < QUALITY; ++i) {
	  std::cout<< i << "%\r";
		shortestDistance = fhDijkstra(seed);          // get the distances of this seed
		for ( int j = 0; j < n ; ++j){                  // associate the vertecies
			agd[j] += shortestDistance[j];      // and calculate agd for each vertex
			if( shortestDistance[j] < mgd[j])
				mgd[j] = shortestDistance[j];
		}
		seed = globalMaxima(mgd);                   // choose new seed vertex		
		delete[] shortestDistance;
	}

	// get the first sample
	int maxInd = globalMaxima(agd);
	
	delete [] agd;
	delete [] mgd;
	std::cout << "100%"<< std::endl;	
	return maxInd;
}
int * HW1::FPS(){
	float * mgd = new float[n];
	int * sample = new int[SAMP_SIZE];
	
	sample[0] = AGD();
	// associate for the first sample and choose the second
	
	mgd = fhDijkstra(sample[0]);
	sample[1] = globalMaxima(mgd);

	float * shortestDistance;
	// we have the first two points
	// sample[0] has all of the associations at first
	// find them for every other sample
	std::cout << std::endl;
	std::cout << "FPS samples: " << std::endl;
	for (int i = 1; i < SAMP_SIZE-1; ++i) {
	  std::cout<< i*100/SAMP_SIZE << "%\r";
		shortestDistance = fhDijkstra(sample[i]); // get the distances of this sample

		for ( int j = 0; j < n ; ++j)               // associate the vertecies
			if(shortestDistance[j] < mgd[j])
				mgd[j] = shortestDistance[j];
		
		sample[i+1] = globalMaxima(mgd);         // choose new sample vertex
		
		delete[] shortestDistance;
	}
	delete[] mgd;
	std::cout << "100%"<< std::endl;
	return sample;
}
int * HW1::sym_inv(){
	/// BASIC FPS
	int * s1_2 = new int[SAMP_SIZE+EXTR_SIZE];
	float * mgd = new float[n];
	s1_2[0] = AGD();
	mgd = fhDijkstra(s1_2[0]);
	s1_2[1] = globalMaxima(mgd);
	float* shortestDistance;
	std::cout << "FPS samples: " << std::endl;
	for (int i = 1; i < SAMP_SIZE; ++i) {
	        s1_2[i] = globalMaxima(mgd);         // choose new sample vertex
		std::cout << i*100/SAMP_SIZE << "%\r";
		shortestDistance = fhDijkstra(s1_2[i]); // get the distances of this sample
		for ( int j = 0; j < n ; ++j)               // associate the vertecies
		 	if( shortestDistance[j] < mgd[j]*0.35)
		 		mgd[j] = shortestDistance[j];
		delete[] shortestDistance;
	}
	///
	std::cout << "100%"<< std::endl;
	int sIdx = SAMP_SIZE;
	int maxMgd = globalMaxima(mgd);
	while (true) {
	  std::cout << "Symmetry-Invariant Samples: " << std::endl;
		for (int i = 0; i < n && sIdx < SAMP_SIZE + EXTR_SIZE; ++i) {
			if(localMaxima(mgd, i) && mgd[i] > mgd[maxMgd]*TOLERANS){
				s1_2[sIdx++] = i;
				std::cout<< i*100/EXTR_SIZE << "%\r";
				shortestDistance = fhDijkstra(s1_2[sIdx-1]);
				for ( int j = 0; j < n ; ++j)    // associate the vertecies
					if( mgd[j] < tolerans(mgd,i)*TOLERANS && shortestDistance[j] < mgd[j])
						mgd[j] = shortestDistance[j];
			}
		}
		if(sIdx < SAMP_SIZE + EXTR_SIZE)
			mgd[maxMgd] *= 0.8;
		else{
			std::cout<< "100%" << std::endl;
			return s1_2;
		}
	}
}
float HW1::tolerans(float * array, int id){
	float sum = 0;
	auto temp = cm->verts[id]->vertList;
	for (auto vert: temp)
		if(array[vert] >= array[id])
			sum+=0;
	return sum/temp.size();
}
bool HW1::localMaxima(float * array, int id){
	for (auto vert: cm->verts[id]->vertList) 
		if(array[vert] >= array[id])
			return false;
	return true;
}
int HW1::globalMaxima(float * funcArray){
	int maxIdx = 0;
	float maxVal = 0;
	int count =0;
	for (int j = 0; j < n ; ++j) {
		if( maxVal < funcArray[j]){
			maxVal = funcArray[j];
			maxIdx = j;
		}
		if(funcArray[j] == 0)
			count++;
	}
	return maxIdx;
}

/// DIJKSTRA IMPLEMENTATIONS
//
float* HW1::aDijkstra(int startingIdx) {

	//      int * prev = new int[n];
        float * dist = new float[n];
        bool * isVisited = new bool[n];

        for (int i = 0; i < n; i++){
		dist[i] = INF;
		//	prev[i] = UNDF;
		isVisited[i] = false;
        }
        dist[startingIdx] = 0;

        int minInd;
        float minValue;
        for (int i = 0; i < n; i++) {

		minValue = INF;
		for (int j = 0; j < n; j++) {
			if (!isVisited[j] && dist[j] < minValue) {
				minInd = j;
				minValue = dist[j];
			}
		}
		isVisited[minInd] = true;
		for (auto neighbor : (cm->verts[minInd]->vertList)){
			float distance = dist[minInd] + cm->getLength(minInd, neighbor); // neighbor Dist
			if(distance < dist[neighbor]){
				dist[neighbor] = distance;
//				prev[neighbor] = minInd;
			}
		}
        }
	return dist;
}
distpair * HW1::hDijkstra(int startingIdx){


	distpair * dist = new distpair[n];
        for (int i = 0; i < n; i++){
		if (i == startingIdx)
			dist[i].y = 0;
		else
			dist[i].y = INF;
		dist[i].x = i;
	}
	
	minHeap<distpair> vertHeap(dist, n);
	
	while(!vertHeap.isEmpty()){
		distpair currentIdx = vertHeap.deleteMin();
		for ( auto neighbor : cm->verts[currentIdx.x]->vertList ){
			float distance = currentIdx.y + cm->getLength(currentIdx.x, neighbor); 
			if( distance < dist[neighbor].y ){
				dist[neighbor].y = distance;
				vertHeap.findAndPercUp(neighbor, distance);
				
			}
		}
	}
	return dist;
}
float* HW1::fhDijkstra(int startingIdx){
	FibonacciHeap<distpair> toExplore;
	node<distpair>** nodePtr = new node<distpair>*[n];
	float* dist = new float[n];

	for (int i =0; i < n; ++i) {
		dist[i] = INF;
		if(i != startingIdx)
			nodePtr[i] = toExplore.insert(distpair(i, INF));
		else
			nodePtr[i] = toExplore.insert(distpair(i, 0));
	}
	dist[startingIdx] = 0;
	
	while(!toExplore.isEmpty()){
		distpair v = toExplore.removeMinimum();
		for (auto neighbor: cm->verts[v.x]->vertList ){
			float distance = dist[v.x] + cm->getLength(v.x, neighbor);
			if( distance < dist[neighbor]){
				dist[neighbor] = distance;
				toExplore.decreaseKey(nodePtr[neighbor], distpair(neighbor, distance));
				
				distpair tmp = nodePtr[neighbor]->getValue();
				auto tmp2 = nodePtr[neighbor];
				nodePtr[neighbor] = nodePtr[tmp.x];
				nodePtr[tmp.x] = tmp2;
			}
		}
	}

	delete[] nodePtr;
	return dist;
}

int* HW1::fhPathDijkstra(int startingIdx, int finishingIdx) {
	
	FibonacciHeap<distpair> toExplore;
	node<distpair>** nodePtr = new node<distpair> * [n];
	float* dist = new float[n];
	int* prevs = new int[n];


	for (int i = 0; i < n; ++i) {
		dist[i] = INF;
		prevs[i] = UNDF;

		if (i != startingIdx)
			nodePtr[i] = toExplore.insert(distpair(i, INF));
		else
			nodePtr[i] = toExplore.insert(distpair(i, 0));
	}
	
	dist[startingIdx] = 0;
	while (!toExplore.isEmpty()) {
		distpair v = toExplore.removeMinimum();
		if (v.x == finishingIdx)
			break;
		for (auto neighbor : cm->verts[v.x]->vertList) {
			float distance = dist[v.x] + cm->getLength(v.x, neighbor);
			if (distance < dist[neighbor]) {
				dist[neighbor] = distance;
				toExplore.decreaseKey(nodePtr[neighbor], distpair(neighbor, distance));
				
				distpair tmp = nodePtr[neighbor]->getValue();
				auto tmp2 = nodePtr[neighbor];
				nodePtr[neighbor] = nodePtr[tmp.x];
				nodePtr[tmp.x] = tmp2;

				prevs[neighbor] = v.x;
			}
		}
	}
	delete[] nodePtr;
	delete[] dist;
	
	return prevs;
}

// WRITE TO FILE
void HW1::write(){
	float** NxN = new float*[n];
	
	for (int i = 0; i < n; ++i) 
		NxN[i] = aDijkstra(i);
	ofstream file;
	file.open("arrayDijMatrix.txt");
	for (int i = 0; i < n; ++i){
		for (int j =0; j<n-1; ++j)
			file << NxN[i][j] << " ";
		file << NxN[i][n-1] << "\n";
	}
	file.close();
	std::cout << "NxN matrix output for the ARRAY DIJKSTRA implementation is written on file." << std::endl;
	distpair ** NxNHeap= new distpair*[n];

	for (int i = 0; i < n; ++i) 
		NxNHeap[i] = hDijkstra(i);
	file.open("heapDijMatrix.txt");
	for (int i = 0; i < n; ++i){
		for (int j =0; j<n-1; ++j)
			file << NxNHeap[i][j].y << " ";
		file << NxNHeap[i][n-1].y << "\n";
	}
	file.close();
	for (int i = 0; i < n; ++i) 
			delete [] NxNHeap[i];
		delete [] NxNHeap;
	std::cout << "NxN matrix output for the HEAP DIJKSTRA implementation is written on file." << std::endl;

	for (int i = 0; i < n; ++i)
		NxN[i] = fhDijkstra(i);
	file.open("fiboHeapDijMatrix.txt");
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n - 1; ++j)
			file << NxN[i][j] << " ";
		file << NxN[i][n - 1] << "\n";
	}
	file.close();
	for (int i = 0; i < n; ++i)
		delete[] NxN[i];
	delete[] NxN;
	std::cout << "NxN matrix output for the FIBONACCI HEAP DIJKSTRA implementation is written on file." << std::endl;
}

// TIMINGS FOR DIJKSTRA IMPLEMENTATIONS
void HW1::timings() {
	std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
	aDijkstra(PATH_START);
	std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	std::cout << "Array Implementation Time: " << time_span.count() << "\n";

	startTime = std::chrono::high_resolution_clock::now();
	hDijkstra(PATH_START);
	endTime = std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	std::cout << "Heap Implementation Time: " << time_span.count() << "\n";

	startTime = std::chrono::high_resolution_clock::now();
	fhDijkstra(PATH_START);
	endTime = std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	std::cout << "Fib Heap Implementation Time: " << time_span.count() << "\n";

}

///   VISUALIZATION   ///
//                        
void HW1::Visualize() {
	// points for sampling q2
	fpsVisualize();
	// points for sampling q3
	symVisualize();
	// edges for dijkstra
	dijVisualize();
	std::cout << "End of calculations.." << std::endl;
}
void HW1::fpsVisualize(){

	std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
	int* samplePoints = FPS();
	std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	std::cout << "Done In " << time_span.count() << "Second" << std::endl;

	
	for (int i = 0; i < SAMP_SIZE; ++i)
		cm->samples.push_back(samplePoints[i]);
	delete [] samplePoints;
	std::cout << "FPS Visualization Complete.." << std::endl;
}
void HW1::dijVisualize(){
	int* path = fhPathDijkstra( PATH_START, PATH_END);
	for (int i = PATH_END; i != PATH_START  ;) {
		for (auto idx:cm->verts[i]->edgeList){
			if (cm->edges[idx]->v1i == path[i] || cm->edges[idx]->v2i == path[i]) {
				cm->sedges.push_back( new Edge(0, i, path[i]) );
				i = path[i];
				break;
			}
		}
	}

	std::cout << "Shortest path visualized using Fibonacci Heap Djikstra implementation.." << std::endl;
	delete[] path;
}
void HW1::symVisualize(){

	std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
	int* samplePoints = sym_inv();
	std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	std::cout << "DONE IN " << time_span.count() << " SECONDS" << std::endl;
	
	for (int i = 0; i < SAMP_SIZE+EXTR_SIZE; ++i) 
		cm->samples.push_back(samplePoints[i]);

	delete [] samplePoints;
}
