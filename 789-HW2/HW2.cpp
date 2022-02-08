#include "HW2.h"
#ifdef HW_2

void HW2::hw2Main(){
#ifdef FINAL

	std::cout << "HOMEWORK ASSIGNMENT 2 " << std::endl << std::endl;
	
	
choice:
	std::cout << "0)  Base Mesh" << std::endl;
	std::cout << "1)  A* Algorithm - Dijkstra  " << std::endl;
	std::cout << "2)  SquareRoot3 - Kobbolt Subdivision " << std::endl;
	std::cout << "3)  4-1 Subdivision " << std::endl;
	std::cout << "4)  Phong Tesselation " << std::endl << std::endl;
	std::cout << "|  info: subdivisions will be applied to \"assets/mesh.off\" " << std::endl;
	std::cout << "|  (q)uit.. " << std::endl << std::endl;
	
	std::cout << "-> ";
	char a;
	std::cin >> a;
	if (a == 'q')
	{
		exit = true;
		goto finish;
	}
	else if (a == '0')
	{
		float sArea = surfArea(cm);
		std::cout << " Base Mesh " << cm->verts.size() << std::endl;
		std::cout << " Vertex: " << cm->verts.size() << std::endl <<
			" Edge: " << cm->edges.size() << std::endl <<
			" Triangle: " << cm->tris.size() << std::endl <<
			" Total Surface Area: " << sArea << std::endl;
		goto finish;
	}
	else if (a == '1')
	{
		astar = true;
	a_ALG:
		std::cout << std::endl;
		std::cout << "    1)  Comparison" << std::endl;
		std::cout << "    2)  Random Distance" << std::endl;
		std::cout << "    3)  A* Distance \\w Input" << std::endl;
		std::cout << "    4)  HW1::Dijkstra \\w Input" << std::endl << std::endl;
		std::cout << "    | (b)ack" << std::endl << std::endl;
		std::cout << "- -> ";
		std::cin >> a;
		if (a == '1') {
			std::cout << std::endl;
			compare();
			goto finish;
		}
		else if (a == '2') {
			srand(time(NULL));
			std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
			aStar(rand() % n, rand() % n);
			std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
			std::cout << "Algorithm applied in " << time_span.count() << std::endl << std::endl;
			goto finish;
		}
		else if (a == '3') {
			int s; int t;
			double len;
  			std::cout << "      StartingIndex: ";
			std::cin >> s;
			std::cout << "      TargetIndex: ";
			std::cin >> t;
			std::cout << std::endl << std::endl;
			srand(time(NULL));
			std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
			aStar(s % n, t % n);
			std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
			len = 0;
			for (int i = 0; i < cm->sedges.size(); ++i)
				len += cm->getLength(cm->sedges[i]->v1i, cm->sedges[i]->v2i); // OR MAYBE PATHLENGTH
			std::cout << "A* Time: " << time_span.count() << "\n";
			std::cout << "Path Length: " << len << std::endl << std::endl;
			goto finish;
		}
		else if (a == '4') {
			int s; int t;
			double len;
			std::cout << "      StartingIndex: ";
			std::cin >> s;
			std::cout << "      TargetIndex: ";
			std::cin >> t;
			std::cout << std::endl << std::endl;
			srand(time(NULL));
			std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
			fhPathDijkstra(s % n, t % n);
			std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
			len = 0;
			for (int i = 0; i < cm->sedges.size(); ++i)
				len += cm->getLength(cm->sedges[i]->v1i, cm->sedges[i]->v2i); // OR MAYBE PATHLENGTH
			std::cout << "Fib Heap Djikstra Time: " << time_span.count() << "\n";
			std::cout << "Path Length: " << len << std::endl << std::endl;
			goto finish;
		}
		else if (a == 'b') {
			goto choice;
		}
		else {
			std::cout << "  Wrong Input try again..." << std::endl << std::endl;
			goto a_ALG;
		}
	}
	else if (a == '2') {
		subdiv = true;
	b_ALG:
		std::cout << std::endl;
		std::cout << "    1)  Just subdivide" << std::endl;
		std::cout << "    2)  Paint vertecies" << std::endl;
		std::cout << "    | (b)ack" << std::endl << std::endl;
		std::cout << "- -> ";
		std::cin >> a;
		if (a == '1') {
		a_subdivide:
			std::cout << std::endl;
			std::cout << "        Subdivision Amount: ";
			int sub;
			std::cin >> sub;
			std::cout << std::endl;
			for (size_t i = 0; i < sub; i++) {
				std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
				nm = _3Subdiv();
				std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
				delete cm;
				cm = nm;
				n = nm->verts.size();
				std::cout << " sqrt(3) Subdivision no: " << i + 1 << " applied in " << time_span.count() << " seconds" << std::endl;
			}


			float sArea = surfArea(nm);
			std::cout << " Vertex: " << nm->verts.size() << std::endl <<
				" Edge: " << nm->edges.size() << std::endl <<
				" Triangle: " << nm->tris.size() << std::endl <<
				" Total Surface Area: " << sArea << std::endl;
			goto finish;
		}
		else if (a == '2') {
			paint = true;
			goto a_subdivide;
		}
		else if (a == 'b') {
			goto choice;
		}
		else {
			std::cout << "  Wrong Input try again..." << std::endl << std::endl;
			goto b_ALG;
		}

		std::cout << std::endl << std::endl;
	}
	else if (a == '3') {
		subdiv = true;
	c_ALG:
		std::cout << std::endl;
		std::cout << "  1)  Just subdivide" << std::endl;
		std::cout << "  2)  Paint vertecies" << std::endl;
		std::cout << "    | (b)ack" << std::endl << std::endl;
		std::cout << "- -> ";
		std::cin >> a;
		if (a == '1') {
		c_subdivide:
			std::cout << std::endl;
			std::cout << "     Subdivision Amount: ";
			int sub;
			std::cin >> sub;
			std::cout << std::endl;
			for (size_t i = 0; i < sub; i++) {
				std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
				nm = _4Subdiv();
				std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
				delete cm;
				cm = nm;
				n = nm->verts.size();
				std::cout << " 4-1 Subdivision no: " << i + 1 << " applied in " << time_span.count() << " seconds" << std::endl;
			}

			float sArea = surfArea(nm);
			std::cout << " Vertex: " << nm->verts.size() << std::endl <<
				" Edge: " << nm->edges.size() << std::endl <<
				" Triangle: " << nm->tris.size() << std::endl <<
				" Total Surface Area: " << sArea << std::endl;
			goto finish;
		}
		else if (a == '2') {
			paint = true;
			goto c_subdivide;
		}
		else if (a == 'b') {
			goto choice;
		}
		else {
			std::cout << "  Wrong Input try again..." << std::endl << std::endl;
			goto c_ALG;
		}

		std::cout << std::endl << std::endl;
	}
	else if (a == '4') {
		subdiv = true;
	d_ALG:
		std::cout << std::endl;
		std::cout << "  1)  Just subdivide" << std::endl;
		std::cout << "  2)  Paint vertecies" << std::endl;
		std::cout << "    | (b)ack" << std::endl << std::endl;
		std::cout << "- -> ";
		std::cin >> a;
		if (a == '1') {
		d_subdivide:
			std::cout << std::endl;
			std::cout << "     Alpha Amount: ";
			std::cin >> alpha;
			std::cout << std::endl;
			std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
			nm = _PSubdiv();
			std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
			delete cm;
			cm = nm;
			n = nm->verts.size();
			std::cout << " Phong Tesselation applied in " << time_span.count() << " seconds" << std::endl;

			float sArea = surfArea(nm);
			std::cout << " Vertex: " << nm->verts.size() << std::endl <<
				" Edge: " << nm->edges.size() << std::endl <<
				" Triangle: " << nm->tris.size() << std::endl <<
				" Total Surface Area: " << sArea << std::endl;
			goto finish;
		}
		else if (a == '2') {
			paint = true;
			goto d_subdivide;
		}
		else if (a == 'b') {
			goto choice;
		}
		else {
			std::cout << "  Wrong Input try again..." << std::endl << std::endl;
			goto d_ALG;
		}

		std::cout << std::endl << std::endl;
	}else{
		std::cout << "Please Choose Another Option.." << std::endl;
		goto choice;
	};
	goto finish;	
	
#endif

#ifdef A_STAR_CROSS 
 a_sc:
	srand(time(NULL));
	aStarCross(rand()%n, rand()%n);
#endif
	
#ifdef A_STAR 
 a_s:	
	srand(time(NULL));
	aStar(rand()%n, rand()%n);
#endif
#ifdef SUB_DIV_3 
 sub_3:
	for (size_t i = 0; i < SUBAMOUNT; i++)
		{
			nm = _3Subdiv();
			delete cm;
			cm = nm;

			n = nm->verts.size();
		}
#endif

#ifdef SUB_DIV_4
 sub_4:
	for (size_t i = 0; i < SUBAMOUNT; i++)
		{
			nm = _4Subdiv();
			delete cm;
			cm = nm;
		
			n = nm->verts.size();
		}
#endif

#ifdef SUB_DIV_P
 sub_p:
	for (size_t i = 0; i < SUBAMOUNT; i++)
		{
			nm = __PSubdiv();
			delete cm;
			cm = nm;
		
			n = nm->verts.size();
		}
#endif

 finish:
	return;
}

double HW2::surfArea(Mesh * mesh){
	double sum = 0; 
	for (int i = 0; i < mesh->tris.size(); ++i) {
		Vec3 v1 = mesh->verts[mesh->tris[i]->v2i]->coords - mesh->verts[mesh->tris[i]->v1i]->coords;
		Vec3 v2 = mesh->verts[mesh->tris[i]->v3i]->coords - mesh->verts[mesh->tris[i]->v1i]->coords;
		sum += (v1^v2).magnitude()/2;
	}
	return sum;

}
int* HW2::fhPathDijkstra(int startingIdx, int finishingIdx) {

	FibonacciHeap<distpair> toExplore;
	node<distpair>** nodePtr = new node<distpair> *[n];
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
		cm->samples.push_back(v.x);
		if (v.x == finishingIdx) {
			while (finishingIdx != startingIdx) {
				cm->sedges.push_back(new Edge(0, finishingIdx, prevs[finishingIdx]));
				finishingIdx = prevs[finishingIdx];
			}

			break;
		}
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
void HW2::compare(){
	srand(time(NULL));
	int k = rand()%n;
	int l = rand()%n;
	double len;

	auto startTime = std::chrono::high_resolution_clock::now();
	fhPathDijkstra(k, l);
	auto endTime = std::chrono::high_resolution_clock::now();
	auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	len = 0;
	for (int i = 0; i < cm->sedges.size() ; ++i) 
		len+= cm->getLength(cm->sedges[i]->v1i, cm->sedges[i]->v2i); // OR MAYBE PATHLENGTH
	std::cout << "Fib Heap Djikstra Time: " << time_span.count() << "\n";
	std::cout << "Path Length: " << len << std::endl<< std::endl;

	cm->sedges.clear();
	cm->samples.clear();
	
	auto startTime1 = std::chrono::high_resolution_clock::now();
	aStar(k, l);
	auto endTime1 = std::chrono::high_resolution_clock::now();
	auto time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(endTime1 - startTime1);
	len = 0;
	for (int i = 0; i <cm->sedges.size() ; ++i) 
		len+=cm->getLength(cm->sedges[i]->v1i, cm->sedges[i]->v2i);
	std::cout << "A* Time: " << time_span1.count() << "\n";
	std::cout << "Path Length: " << len << std::endl<< std::endl;
	
}
double HW2::pathLength(int* path, int s, int t ){
	double sum = 0;
	for (int i = t; i != s  ;) {
		for (auto idx:cm->verts[i]->edgeList){
			if (cm->edges[idx]->v1i == path[i] || cm->edges[idx]->v2i == path[i]) {
				sum += cm->getLength(i, path[i]);
				i = path[i];
				break;
			}
		}
	}
	delete[] path;
	return sum;
}



void HW2::aStarCross(int s, int t) {
	float* gScore = new float[n];
	float* fScore = new float[n];
	bool* isInside = new bool[n];

	int* cameFrom = new int[n];

	for (int i = 0; i < n; i++) {
		gScore[i] = INF;
		fScore[i] = INF;
		isInside[i] = false;
	}
	gScore[s] = 0;
	fScore[s] = h(s, t);
	isInside[s] = true;

	distpair* pairs = new distpair[1];
	pairs[0] = distpair(s, fScore[s]);
	minHeap<distpair> vertHeap(pairs, 1);

	while (!vertHeap.isEmpty()) {
		distpair currentIdx = vertHeap.deleteMin();
		cm->samples.push_back(currentIdx.x);
		isInside[currentIdx.x] = false;
		if (currentIdx.x == t) {
			delete[] gScore;
			delete[] fScore;
			delete[] isInside;
			delete[] pairs;

			while (t != s) {
				cm->sedges.push_back(new Edge(0, t, cameFrom[t]));
				t = cameFrom[t];
			}
			delete[] cameFrom;
			return;
		}

		for (auto neighbor : cm->verts[currentIdx.x]->vertList) {
			float t_g = gScore[currentIdx.x] + cm->getLength(currentIdx.x, neighbor);
			if (t_g < gScore[neighbor]) {
				cameFrom[neighbor] = currentIdx.x;
				gScore[neighbor] = t_g;
				fScore[neighbor] = gScore[neighbor] + h(neighbor, t);
				if (!isInside[neighbor]) {
					distpair pair(neighbor, fScore[neighbor]);
					vertHeap.insert(pair);
					isInside[neighbor] = true;
				}

			}
		}
	}
	delete[] gScore;
	delete[] fScore;
	delete[] isInside;
	delete[] cameFrom;
	delete[] pairs;

	throw "COULD NOT FIND THE TARGET NODE";
}

void HW2::aStar(int s, int t){
	float * gScore = new float[n];
	float * fScore = new float[n];
	bool * isInside = new bool[n];
	
	int * cameFrom = new int[n];

	for (int i = 0; i < n; i++){
		gScore[i] = INF;
		fScore[i] = INF;
		isInside[i] = false;
	}
	gScore[s] = 0;
	fScore[s] = h(s, t);
	isInside[s] = true;
	
	distpair* pairs = new distpair[1];
	pairs[0] = distpair(s, fScore[s]);
	minHeap<distpair> vertHeap(pairs, 1);
	
	while(!vertHeap.isEmpty()){
		distpair currentIdx = vertHeap.deleteMin();
		cm->samples.push_back(currentIdx.x);
		isInside[currentIdx.x] =  false;
		if (currentIdx.x == t) {
			delete [] gScore;
			delete [] fScore;
			delete [] isInside;
			delete [] pairs;

			while(t != s){
				cm->sedges.push_back( new Edge(0, t, cameFrom[t]));
				t = cameFrom[t];
			}
			delete [] cameFrom;
			return;
		}
		
		for ( auto neighbor : cm->verts[currentIdx.x]->vertList ){
			float t_g = gScore[currentIdx.x] + cm->getLength(currentIdx.x, neighbor); 
			if( t_g < gScore[neighbor] ){
				cameFrom[neighbor] = currentIdx.x;
				gScore[neighbor] = t_g;
				fScore[neighbor] = gScore[neighbor] + h(neighbor,t);
				if(!isInside[neighbor]){
					distpair pair(neighbor, fScore[neighbor]);
					vertHeap.insert(pair);
					isInside[neighbor] = true;
				}
				
			}
		}
	}
	delete [] gScore;
	delete [] fScore;
	delete [] isInside;
	delete [] cameFrom;
	delete [] pairs;
	
	throw "COULDNT FIND THE TARGET NODE";
}

float HW2::h(int s, int t){
	return (cm->getLength(s, t));
}


////////
//// Q2

Mesh* HW2::_3Subdiv() {
	Mesh* sub = new Mesh;
	int vertSize = cm->verts.size();

	for (int i = 0; i < vertSize; ++i) {
		Vertex* vert = new Vertex(i, cm->verts[i]);
		sub->verts.push_back(vert);
	}

	// ADD VERTICIES
	vector<float> dists;
	vector<int> colIds;
	int ii;
	float mindist = INF;
	for (auto tri : cm->tris) {
		Vec3 midpoint = (cm->verts[tri->v1i]->coords +
			cm->verts[tri->v2i]->coords +
			cm->verts[tri->v3i]->coords) / 3;

		sub->addVertex(midpoint.x, midpoint.y, midpoint.z);
		ii = sub->verts.size() - 1;
		colIds.push_back(ii);
		float d = (sub->getLength(tri->v1i, ii) + sub->getLength(tri->v2i, ii) + sub->getLength(tri->v3i, ii));
		if (mindist > d)
			mindist = d;
		dists.push_back(d);

	}
	float maxdist = -1;
	for (int i = 0; i < dists.size(); ++i) {
		dists[i] = dists[i] / mindist - 1;
		if (maxdist < dists[i])
			maxdist = dists[i];
	}
	for (int i = 0; i < dists.size(); ++i) {
		float a = (dists[i] / maxdist) * 255;
		sub->verts[colIds[i]]->color = SbVec3f(a + 0, 255 - a, 0);
	}

	// SMOOTHING
	for (int i = 0; i < vertSize; ++i) {
		int n = cm->verts[i]->vertList.size();
		float an = (4 - 2*(cos(2*M_PI/n)))/9;

		Vec3 sum(0,0,0);
		for (auto neigh : cm->verts[i]->vertList)
			sum = sum + cm->verts[neigh]->coords;

		sub->verts[i]->coords = (sub->verts[i]->coords*(1-an)) + (sum*(an/n));
	}
	
	// CONNECT & FLIP
	int edgeSize = cm->edges.size();
	bool* isHandl = new bool[edgeSize];
	for (int ed = 0; ed < edgeSize; ed++)
		{
			isHandl[ed] = false;
		}


	int faceSize = cm->tris.size();
	for (int i = 0; i < faceSize; ++i) {
		vector<int> currentEdges = getEdges(cm->tris[i]);
		for (int j = 0; j < 3; j++) {
			int ed = currentEdges[j];
			//cm->sedges.push_back(cm->edges[ed]);

			int sing = cm->edges[ed]->isNotSingle();
			if (sing) {              // 1 or tris
				isHandl[ed] = true;
				if (sing == 1) {                //sagda
					sub->addTriangle( 	  // v1 sag v2
							 (cm->edges[ed]->v1i),
							 (cm->edges[ed]->t1i + vertSize),
							 (cm->edges[ed]->v2i));
				}
				else if (sing == 2) {               //solda
					sub->addTriangle(         // v1 v2 sol
							 (cm->edges[ed]->v1i),
							 (cm->edges[ed]->v2i),
							 (cm->edges[ed]->t2i + vertSize));
				}
			}
			else if(!isHandl[ed]) { // 2 tris
				isHandl[ed] = true;
				// sol sag ileri
				sub->addTriangle(
						 (cm->edges[ed]->t2i + vertSize),
						 (cm->edges[ed]->t1i + vertSize),
						 (cm->edges[ed]->v2i));
				//  geri sag sol
				sub->addTriangle(
						 (cm->edges[ed]->v1i),
						 (cm->edges[ed]->t1i + vertSize),
						 (cm->edges[ed]->t2i + vertSize));
			}
		}
	}
	return sub;
}



Mesh* HW2::_4Subdiv() {

	// CREATE
	Mesh* sub = new Mesh;
	int vertSize = cm->verts.size();

	for (int i = 0; i < vertSize; ++i) {
		Vertex* vert = new Vertex(i, cm->verts[i]);
		sub->verts.push_back(vert);
	}

	// ADD VERTICIES
	vector<float> dists;
	vector<int> colIds;
	int ii;
	int edgeSize = cm->edges.size();
	float mindist = INF;
	for (int ed = 0; ed < edgeSize; ++ed) {
		sub->addVertex((cm->verts[cm->edges[ed]->v1i]->coords + cm->verts[cm->edges[ed]->v2i]->coords) / 2);

		ii = sub->verts.size() - 1;
		colIds.push_back(ii);
		float d = sub->getLength(cm->edges[ed]->v1i, cm->edges[ed]->v2i);
		if (mindist > d)
			mindist = d;
		dists.push_back(d);
	}
	float maxdist = -1;
	for (int i = 0; i < dists.size(); ++i) {
		dists[i] = dists[i] / mindist - 1;
		if (maxdist < dists[i])
			maxdist = dists[i];
	}
	for (int i = 0; i < dists.size(); ++i) {
		float a = (dists[i] / maxdist) * 255;
		sub->verts[colIds[i]]->color = SbVec3f( a + 0, 255-a, 0);
	}
	
	// ADD TRIANGLES (CONNECT & FLIP)
	int faceSize = cm->tris.size();
	for (int i = 0; i < faceSize; ++i) {
		vector<int> currentEdges = getEdges(cm->tris[i]);		
		int v1 = cm->tris[i]->v1i;
		int v2 = cm->tris[i]->v2i;
		int v3 = cm->tris[i]->v3i;
		
		// v1 e1 e3
		sub->addTriangle( v1,
				  (currentEdges[0] + vertSize),
				  (currentEdges[2] + vertSize));

		// e1 v2 e2
		sub->addTriangle( (currentEdges[0] + vertSize),
				  v2,
				  (currentEdges[1] + vertSize));

		// e3 e1 e2a
		sub->addTriangle( (currentEdges[2] + vertSize),
				  (currentEdges[0] + vertSize),
				  (currentEdges[1] + vertSize));
		
		// e3 e2 v3
		sub->addTriangle( (currentEdges[2] + vertSize),
				  (currentEdges[1] + vertSize),
				  v3);
	}
	return sub;
}

Vec3 HW2::tangentPlaneCorrespondence(Vec3 p, Vec3 n, Vec3 c){
	return  Vec3(n.x*(c.x-p.x), n.y*(c.y-p.y), n.z*(c.z-p.z));
}

vector<Vec3> HW2::createPointSet(Triangle* tri, int beta, vector<int>pMods){ // vector<Vertex*>
	vector<Vec3> bary;
	Vec3 ori = cm->coord(tri->v1i);
	Vec3 v1  = cm->coord(tri->v2i)-ori;
	Vec3 v2  = cm->coord(tri->v3i)-ori;
	float v1m = v1.magnitude();
	float v2m = v2.magnitude();

	// bary.push_back(ori); // for the  barycentric
	// bary.push_back(v1);  // navigation    inside
	// bary.push_back(v2);  // the current triangle

	// bary.push_back(Vec3(0,0,1)); // v1i
	// bary.push_back(Vec3(1,0,0)); // v2i
	// bary.push_back(Vec3(0,1,0)); // v3i
	

	int curEdgeMod;
	int div;
	int edgeNo;
	vector<int> ed = getEdges(tri);
	for (int i = 0; i < 3; ++i) {
		Vec3 start;
		Vec3 divesee;
		if( (cm->edges[ed[i]]->v1i == tri->v2i && cm->edges[ed[i]]->v2i == tri->v1i) ||
		    (cm->edges[ed[i]]->v2i == tri->v2i && cm->edges[ed[i]]->v1i == tri->v1i)   ){ // v1/v2 edge
			start = cm->coord(tri->v1i);
			divesee = cm->coord(tri->v2i) - start;
		}
		else if( (cm->edges[ed[i]]->v1i == tri->v3i && cm->edges[ed[i]]->v2i == tri->v1i) ||
			 (cm->edges[ed[i]]->v2i == tri->v3i && cm->edges[ed[i]]->v1i == tri->v1i)       ){ // v1/v3 edge
			start = cm->coord(tri->v1i);
			divesee = cm->coord(tri->v3i) - start;
		}
		else if( (cm->edges[ed[i]]->v1i == tri->v2i && cm->edges[ed[i]]->v2i == tri->v3i) ||
			 (cm->edges[ed[i]]->v2i == tri->v2i && cm->edges[ed[i]]->v1i == tri->v3i)){ // v2/v3 edge
			start = cm->coord(tri->v2i);
			divesee = cm->coord(tri->v3i) - start;
		}

		
		if(pMods[cm->edges[ed[i]]->t1i] > pMods[cm->edges[ed[i]]->t2i])
			curEdgeMod = pMods[cm->edges[ed[i]]->t1i];
		else
			curEdgeMod = pMods[cm->edges[ed[i]]->t2i];

		switch(curEdgeMod){
		case 1:
			div = 2;
			break;
		case 2:
			div = 2;
			break;
		case 3:
			div = 3;
			break;
		case 4:
			div = 4;
			break;
		default:
			throw "REALLY ???";
		}
		Vec3 bar = (divesee/div);
		for (int j = 1; j < div ; ++j) {
			Vec3 pos = (start + bar*j)-ori;
			Vec3 tmp = Vec3( ((pos*v1)/pow(v1m,2)),((pos*v2)/pow(v2m,2)), 0 );
			tmp.z = 1-(tmp.x + tmp.y);
			bary.push_back(tmp);
		}

		
	}
	
	switch(beta){
	case 1:
		break;
	case 2:
		bary.push_back(Vec3(1.0/3, 1.0/3, 1.0/3));
		break;
	case 3:
		bary.push_back(Vec3((2.0/9), (2.0/9), (5.0/9)));
		bary.push_back(Vec3((5.0/9), (2.0/9), (2.0/9)));
		bary.push_back(Vec3((2.0/9), (5.0/9), (2.0/9)));			
		break;
	case 4:
		bary.push_back(Vec3((1.0/6), (1.0/6), (4.0/6)));
		bary.push_back(Vec3((4.0/6), (1.0/6), (1.0/6)));
		bary.push_back(Vec3((1.0/6), (4.0/6), (1.0/6)));

		bary.push_back(Vec3((5.0/12), (2.0/12), (5.0/12)));
		bary.push_back(Vec3((2.0/12), (5.0/12), (5.0/12)));
		bary.push_back(Vec3((5.0/12), (5.0/12), (2.0/12)));

		bary.push_back(Vec3((1.0/3), (1.0/3), (1.0/3)));
		break;
	default:
		throw "REALLY ???";
	}
	return bary;
}

Mesh * HW2::__PSubdiv(){
	Mesh * sub = _4Subdiv();

	// for (int i = 0 ; i < cm->verts.size(); ++i) {
	// 	// for old verts
	// }
	
	// for (int i = cm->verts.size() ; i < sub->verts.size(); ++i) {
	// 	// for new verts
	// }

	// FIND VERTEX NORMALS
	for (int i = 0 ; i < cm->verts.size(); ++i) {
		cm->verts[i]->normal = Vec3(0,0,0);
		for (int j = 0 ; j < cm->verts[i]->triList.size(); ++j){
			Triangle * tri = cm->tris[cm->verts[i]->triList[j]];
			cm->verts[i]->normal = cm->verts[i]->normal + findNormal(tri->v1i, tri->v2i, tri->v3i);// * cm->getArea(tri);
		}
		cm->verts[i]->normal = cm->verts[i]->normal.normalize(); // normal is set
	}	
	for (int i = 0 ; i < cm->tris.size(); ++i) {
		Triangle* tri = cm->tris[i];
		
		Vec3 v1 = cm->verts[cm->tris[i]->v1i]->coords;
		Vec3 v2 = cm->verts[cm->tris[i]->v2i]->coords;
		Vec3 v3 = cm->verts[cm->tris[i]->v3i]->coords;

		// for old tris
		auto eds = getEdges(tri);
		for (int j = 0; j < 3; ++j) {
			Vertex* point = sub->verts[cm->verts.size() + cm->edges[eds[j]]->idx];
			Vec3 p = point->coords;

			float area = cm->getArea(tri);
			//float w1 = ((v3-p)^(v2-p)).magnitude()/(2* area);
			//float w2 = ((v3-p)^(v1-p)).magnitude()/(2* area);
			//float w3 = ((v2-p)^(v1-p)).magnitude()/(2* area);
			
			float w1 = ((p * v1) / pow(v1.magnitude(), 2));
			float w2 = ((p * v2) / pow(v2.magnitude(), 2)); 
			float w3 = ((p * v3) / pow(v3.magnitude(), 2)); 

			Vec3 bary = Vec3(w1, w2, w3).normalize();
			point->coords = p_star(bary, tri);
		}
	}
	return sub;
}


Mesh * HW2::_PSubdiv(){
	// CREATE
	Mesh * sub = new Mesh;

	int vertSize = cm->verts.size();
	for (int i = 0 ; i < vertSize; ++i) 
		sub->addVertex(cm->verts[i]->coords);

	// FIND VERTEX NORMALS
	for (int i = 0 ; i < cm->verts.size(); ++i) {
		cm->verts[i]->normal = Vec3(0,0,0);
		for (int j = 0 ; j < cm->verts[i]->triList.size(); ++j){
			Triangle * tri = cm->tris[cm->verts[i]->triList[j]];
			cm->verts[i]->normal = cm->verts[i]->normal +
				findNormal(tri->v1i, tri->v2i, tri->v3i) * cm->getArea(tri);
		}
		cm->verts[i]->normal = cm->verts[i]->normal.normalize(); // normal is set
	}

	vector<float> dists;
	vector<int> colIds;
	int ii;
	float mindist = INF;
	for (int i = 0; i < cm->tris.size(); ++i) {

		//sub->addVertex(chooseBary(cm->tris[i]));

		Vec3 uvw = Vec3(1. / 3, 1. / 3, 1. / 3);
		sub->addVertex(p_star(uvw, cm->tris[i]));

		sub->addTriangle(cm->tris[i]->v1i, cm->tris[i]->v2i, sub->verts.size() - 1);
		sub->addTriangle(cm->tris[i]->v2i, cm->tris[i]->v3i, sub->verts.size() - 1);
		sub->addTriangle(cm->tris[i]->v3i, cm->tris[i]->v1i, sub->verts.size() - 1);

		ii = sub->verts.size() - 1;
		colIds.push_back(ii);
		float d = (sub->getLength(cm->tris[i]->v1i, ii) + sub->getLength(cm->tris[i]->v2i, ii) + sub->getLength(cm->tris[i]->v3i, ii));
		if (mindist > d)
			mindist = d;
		dists.push_back(d);

	}
	float maxdist = -1;
	for (int i = 0; i < dists.size(); ++i) {
		dists[i] = dists[i] / mindist - 1;
		if (maxdist < dists[i])
			maxdist = dists[i];
	}
	for (int i = 0; i < dists.size(); ++i) {
		float a = (dists[i] / maxdist) * 255;
		sub->verts[colIds[i]]->color = SbVec3f(a + 0, 255 - a, 0);
	}
	
	return sub;
}

Vec3 HW2::chooseBary(Triangle* tri, float div) {
	Vec3 bar;
	Vec3 tmp;
	Vec3 keep;
	int scan = 1 / div;
	float max = -INF;
	for (size_t i = 0; i < scan; i++)
	{
		for (size_t j = 0; j < scan; j++)
		{
			for (size_t k = 0; k < scan; k++)
			{
				tmp = p_star(bar, tri);
				if (tmp.magnitude() > max) {
					keep = tmp;
					max = tmp.magnitude();
				}
				bar = bar + Vec3(0, 0, div);
			}
			bar = bar + Vec3(0, div, 0);
		}
		bar = bar + Vec3(div, 0, 0);

	}

	return keep;
}

/////// EARLIER PHONG
//Mesh* HW2::_PSubdiv() {
//	// CREATE
//	Mesh* sub = new Mesh;
//
//	int vertSize = cm->verts.size();
//	for (int i = 0; i < vertSize; ++i)
//		sub->addVertex(cm->verts[i]->coords);
//
//	// FIND VERTEX NORMALS
//	for (int i = 0; i < cm->verts.size(); ++i) {
//		cm->verts[i]->normal = Vec3(0, 0, 0);
//		for (int j = 0; j < cm->verts[i]->triList.size(); ++j) {
//			Triangle* tri = cm->tris[cm->verts[i]->triList[j]];
//			cm->verts[i]->normal = cm->verts[i]->normal +
//				findNormal(tri->v1i, tri->v2i, tri->v3i) * cm->getArea(tri);
//		}
//		cm->verts[i]->normal = cm->verts[i]->normal.normalize(); // normal is set
//	}
//
//	////// POINT MOD DECISION
//	// vector<float> areas;
//	// float minArea = INF, maxArea = -1;
//	// int curMin = -1, curMax = -1;
//	// for (int i = 0 ; i < cm->tris.size(); ++i) {
//	// 	areas.push_back(cm->getArea(cm->tris[i]));
//	// 	if (areas[i] < minArea){
//	// 		minArea = areas[i];
//	// 		curMin = i;
//	// 	}
//	// }
//
//	// maxArea = -1;
//	// curMax = -1;
//	// for (int i = 0 ; i < areas.size(); ++i) {
//	// 	areas[i] = (areas[i] / minArea) - 1;
//	// 	if (areas[i] > maxArea){
//	// 		maxArea = areas[i];
//	// 		curMax = i;			
//	// 	}
//	// }
//
//	// vector<int>pointMod;
//	// for (int i = 0 ; i < areas.size(); ++i) {
//	// 	areas[i] = (areas[i] / maxArea)*(MODCOUNT-0.1) + 1;
//	// 	pointMod.push_back((int)areas[i]); // floor
//	// }
//	vector<Vec3> nPointSetBary = basicPointSet(MODCOUNT);
//	for (int i = 0; i < cm->tris.size(); ++i) {
//		vector<Vec3> nPointSetCartesian;
//		Vec3 ori = sub->coord(cm->tris[i]->v1i);
//		Vec3 v1 = sub->coord(cm->tris[i]->v2i) - ori;
//		Vec3 v2 = sub->coord(cm->tris[i]->v3i) - ori;
//
//		vector<Vertex*> npoinn;
//		for (int j = 0; j < nPointSetBary.size(); ++j) {
//			nPointSetCartesian.push_back(ori + v1 * nPointSetBary[j].x + v2 * nPointSetBary[j].y);
//
//			sub->addVertex(nPointSetCartesian[j]);
//			npoinn.push_back(sub->verts[sub->verts.size() - 1]);
//		}
//		npoinn.insert(npoinn.begin() + 0, sub->verts[cm->tris[i]->v3i]);
//		npoinn.insert(npoinn.begin() + 0, sub->verts[cm->tris[i]->v2i]);
//		npoinn.insert(npoinn.begin() + 0, sub->verts[cm->tris[i]->v1i]);
//
//		basicTriangulate(sub, npoinn, MODCOUNT);
//
//
//		for (int j = 3; j < npoinn.size(); j++)
//			npoinn[j]->coords = nPointSetCartesian[j - 3];// p_star(nPointSetBary[j-3], cm->tris[i] );
//
//		nPointSetCartesian.clear();
//	}
//	return sub;
//	//////// CREATE USING COMPLEX METHOD
//	// vector<Vec3> nPointSetBary = createPointSet( cm->tris[i], pointMod[i], pointMod);
//	// vector<Vec3> nPointSetCartesian;
//
//	// vector<Vertex*> npoinn;
//	// for (int j = 0; j < nPointSetBary.size(); ++j) {
//	// 	nPointSetCartesian.push_back( ori + v1*nPointSetBary[j].x + v2*nPointSetBary[j].y);
//
//	// 	sub->addVertex(nPointSetCartesian[j]);
//	// 	npoinn.push_back(sub->verts[sub->verts.size()-1]);
//	// }
//	// npoinn.insert(npoinn.begin() + 0, sub->verts[cm->tris[i]->v3i]);
//	// npoinn.insert(npoinn.begin() + 0, sub->verts[cm->tris[i]->v2i]);
//	// npoinn.insert(npoinn.begin() + 0, sub->verts[cm->tris[i]->v1i]);
//
//	// triangulate(sub, npoinn);
//
//	// for (int j = 3; j < npoinn.size(); j++) 
//	// 	npoinn[j]->coords = p_star(nPointSetBary[j-3], cm->tris[i] );
//
//	// vec3 uvw = Vec3(1./3, 1./3, 1. - 2./3 );
//	// sub->addVertex( p_star( uvw, cm->tris[i] ) );
//
//	// sub->addTriangle( cm->tris[i]->v1i, cm->tris[i]->v2i, sub->verts.size()-1 );
//	// sub->addTriangle( cm->tris[i]->v2i, cm->tris[i]->v3i, sub->verts.size()-1 );
//	// sub->addTriangle( cm->tris[i]->v3i, cm->tris[i]->v1i, sub->verts.size()-1 );
//
//}
void HW2::basicTriangulate(Mesh* sub, vector<Vertex*> pset, int mod){

	switch(mod){
	case 1:
	    sub->addTriangle(pset[3]->idx, pset[0]->idx, pset[5]->idx); // sol v1 alt
		sub->addTriangle(pset[4]->idx, pset[1]->idx, pset[5]->idx); // alt v2 sag >>
		sub->addTriangle(pset[4]->idx, pset[5]->idx, pset[3]->idx); // sol alt sag >>
		sub->addTriangle(pset[3]->idx, pset[4]->idx, pset[2]->idx); // sol sag v3

		/////////////////////////////////////////////////// v1 0
		/////////////////////////////////////////////////// v2 1
		/////////////////////////////////////////////////// v3 2		
		
		// // // // bary.push_back(Vec3(1.0/2, 0, 1.0/2)); // sol 3
		// // // // bary.push_back(Vec3(1.0/2, 1.0/2, 0)); // sag 4
		// // // // bary.push_back(Vec3(0, 1.0/2, 1.0/2)); // alt 5
		break;
	case 2:
		sub->addTriangle(pset[0]->idx, pset[5]->idx, pset[6]->idx);
		sub->addTriangle(pset[5]->idx, pset[1]->idx, pset[6]->idx);

		sub->addTriangle(pset[0]->idx, pset[6]->idx, pset[3]->idx);
		sub->addTriangle(pset[6]->idx, pset[1]->idx, pset[4]->idx);

		sub->addTriangle(pset[3]->idx, pset[6]->idx, pset[2]->idx);
		sub->addTriangle(pset[6]->idx, pset[4]->idx, pset[2]->idx);

		// bary.push_back(Vec3(1.0/2, 0, 1.0/2)); // 3 sol 
		// bary.push_back(Vec3(1.0/2, 1.0/2, 0)); // 4 sag
		// bary.push_back(Vec3(0, 1.0/2, 1.0/2)); // 5 alt

		// bary.push_back(Vec3(1.0/3, 1.0/3, 1.0/3)); // 6 middle
		break;
	case 3:
		sub->addTriangle(pset[0]->idx, pset[5]->idx, pset[9]->idx);
		sub->addTriangle(pset[0]->idx, pset[9]->idx, pset[3]->idx);
		
		sub->addTriangle(pset[5]->idx, pset[10]->idx, pset[9]->idx);
		sub->addTriangle(pset[5]->idx, pset[6]->idx, pset[10]->idx);

		sub->addTriangle(pset[6]->idx, pset[1]->idx, pset[10]->idx);
		sub->addTriangle(pset[1]->idx, pset[7]->idx, pset[10]->idx);

		sub->addTriangle(pset[3]->idx, pset[9]->idx, pset[4]->idx);		
		sub->addTriangle(pset[9]->idx, pset[11]->idx, pset[4]->idx);

		sub->addTriangle(pset[10]->idx, pset[8]->idx, pset[7]->idx);		
		sub->addTriangle(pset[10]->idx, pset[7]->idx, pset[11]->idx);

		sub->addTriangle(pset[4]->idx, pset[11]->idx, pset[2]->idx);
		sub->addTriangle(pset[11]->idx, pset[7]->idx, pset[2]->idx);

		sub->addTriangle(pset[9]->idx, pset[10]->idx, pset[11]->idx);
		// bary.push_back(Vec3((1.0/3),     (0), (2.0/3))); sol1 3
		// bary.push_back(Vec3((2.0/3),     (0), (1.0/3))); sol2 4

		// bary.push_back(Vec3(    (0), (1.0/3), (2.0/3))); alt1 5
		// bary.push_back(Vec3(    (0), (2.0/3), (1.0/3))); alt2 6
		
		// bary.push_back(Vec3((2.0/3), (1.0/3),     (0))); sag2 7
		// bary.push_back(Vec3((1.0/3), (2.0/3),     (0))); sag1 8

		// bary.push_back(Vec3((2.0/9), (2.0/9), (5.0/9))); m1 9
		// bary.push_back(Vec3((5.0/9), (2.0/9), (2.0/9))); m2 10
		// bary.push_back(Vec3((2.0/9), (5.0/9), (2.0/9))); m3 11

		break;
	case 4:
		//row 1
		sub->addTriangle(pset[0]->idx, pset[6]->idx, pset[12]->idx);
		sub->addTriangle(pset[0]->idx, pset[12]->idx, pset[3]->idx);

		sub->addTriangle(pset[6]->idx, pset[7]->idx, pset[16]->idx);// >>>
		sub->addTriangle(pset[6]->idx, pset[16]->idx, pset[12]->idx);
		
		sub->addTriangle(pset[7]->idx, pset[8]->idx, pset[13]->idx);
		sub->addTriangle(pset[7]->idx, pset[14]->idx, pset[16]->idx);

		sub->addTriangle(pset[8]->idx, pset[1]->idx, pset[13]->idx);		
		sub->addTriangle(pset[1]->idx, pset[10]->idx, pset[13]->idx);

		//row 2
		sub->addTriangle(pset[3]->idx, pset[12]->idx, pset[15]->idx);
		sub->addTriangle(pset[3]->idx, pset[15]->idx, pset[4]->idx);

		sub->addTriangle(pset[12]->idx, pset[16]->idx, pset[18]->idx);
		sub->addTriangle(pset[12]->idx, pset[18]->idx, pset[15]->idx);

		sub->addTriangle(pset[16]->idx, pset[13]->idx, pset[18]->idx);
		sub->addTriangle(pset[13]->idx, pset[17]->idx, pset[18]->idx);
		
		sub->addTriangle(pset[13]->idx, pset[9]->idx, pset[10]->idx);
		sub->addTriangle(pset[13]->idx, pset[10]->idx, pset[17]->idx);

		//row 3
		sub->addTriangle(pset[4]->idx, pset[15]->idx, pset[14]->idx);
		sub->addTriangle(pset[4]->idx, pset[14]->idx, pset[5 ]->idx);

		sub->addTriangle(pset[15]->idx, pset[18]->idx, pset[14]->idx);
		sub->addTriangle(pset[18]->idx, pset[17]->idx, pset[14]->idx);

		sub->addTriangle(pset[17]->idx, pset[10 ]->idx, pset[11]->idx);
		sub->addTriangle(pset[17]->idx, pset[11]->idx, pset[14]->idx);

		//row 4
		sub->addTriangle(pset[5 ]->idx, pset[14]->idx, pset[ 2]->idx);
		sub->addTriangle(pset[14]->idx, pset[11]->idx, pset[ 2]->idx);

		/////////////////////////////////////////////////// v1 0
		/////////////////////////////////////////////////// v2 1
		/////////////////////////////////////////////////// v3 2		

		// bary.push_back(Vec4((1.0/4),     (0), (3.0/4))); sol1 3
		// bary.push_back(Vec4((2.0/4),     (0), (2.0/4))); sol2 4
		// bary.push_back(Vec4((3.0/4),     (0), (1.0/4))); sol3 5

		// bary.push_back(Vec4(    (0), (1.0/4), (3.0/4))); alt1 6
		// bary.push_back(Vec4(    (0), (2.0/4), (2.0/4))); alt2 7
		// bary.push_back(Vec4(    (0), (3.0/4), (1.0/4))); alt3 8
		
		// bary.push_back(Vec4((1.0/4), (3.0/4),     (0))); sag1 9
		// bary.push_back(Vec4((2.0/4), (2.0/4),     (0))); sag2 10
		// bary.push_back(Vec4((3.0/4), (1.0/4),     (0))); sag3 11

		// // INSIDE
		// bary.push_back(Vec3((1.0/6), (1.0/6), (4.0/6))); m1 12
		// bary.push_back(Vec3((4.0/6), (1.0/6), (1.0/6))); m2 13
		// bary.push_back(Vec3((1.0/6), (4.0/6), (1.0/6))); m3 14

		// bary.push_back(Vec3((5.0/12), (2.0/12), (5.0/12))); msol1 15
		// bary.push_back(Vec3((2.0/12), (5.0/12), (5.0/12))); malt1 16
		// bary.push_back(Vec3((5.0/12), (5.0/12), (2.0/12))); msag1 17

		// bary.push_back(Vec3((1.0/3), (1.0/3), (1.0/3))); mid 18
		break;
	default:
		throw "REALLY ???";
	};

	
}

vector<Vec3>  HW2::basicPointSet(int mod){
	vector<Vec3> bary;

	switch (mod) {
	case 1:
		// ON EDGES
		bary.push_back(Vec3(1.0 / 2, 0, 1.0 / 2));
		bary.push_back(Vec3(1.0 / 2, 1.0 / 2, 0));
		bary.push_back(Vec3(0, 1.0 / 2, 1.0 / 2));
		break;
	case 2:
		// ON EDGES
		bary.push_back(Vec3(1.0 / 2, 0, 1.0 / 2));
		bary.push_back(Vec3(1.0 / 2, 1.0 / 2, 0));
		bary.push_back(Vec3(0, 1.0 / 2, 1.0 / 2));

		// INSIDE
		bary.push_back(Vec3(1.0 / 3, 1.0 / 3, 1.0 / 3));
		break;
	case 3:
		// ON EDGES
		bary.push_back(Vec3((1.0 / 3), (0), (2.0 / 3)));
		bary.push_back(Vec3((2.0 / 3), (0), (1.0 / 3)));

		bary.push_back(Vec3((0), (1.0 / 3), (2.0 / 3)));
		bary.push_back(Vec3((0), (2.0 / 3), (1.0 / 3)));

		bary.push_back(Vec3((2.0 / 3), (1.0 / 3), (0)));
		bary.push_back(Vec3((1.0 / 3), (2.0 / 3), (0)));

		// INSIDE
		bary.push_back(Vec3((2.0 / 9), (2.0 / 9), (5.0 / 9)));
		bary.push_back(Vec3((5.0 / 9), (2.0 / 9), (2.0 / 9)));
		bary.push_back(Vec3((2.0 / 9), (5.0 / 9), (2.0 / 9)));
		break;
	case 4:

		bary.push_back(Vec3((1.0 / 4), (0), (3.0 / 4)));
		bary.push_back(Vec3((2.0 / 4), (0), (2.0 / 4)));
		bary.push_back(Vec3((3.0 / 4), (0), (1.0 / 4)));

		bary.push_back(Vec3((0), (1.0 / 4), (3.0 / 4)));
		bary.push_back(Vec3((0), (2.0 / 4), (2.0 / 4)));
		bary.push_back(Vec3((0), (3.0 / 4), (1.0 / 4)));

		bary.push_back(Vec3((1.0 / 4), (3.0 / 4), (0)));
		bary.push_back(Vec3((2.0 / 4), (2.0 / 4), (0)));
		bary.push_back(Vec3((3.0 / 4), (1.0 / 4), (0)));

		// INSIDE
		bary.push_back(Vec3((1.0 / 6), (1.0 / 6), (4.0 / 6)));
		bary.push_back(Vec3((4.0 / 6), (1.0 / 6), (1.0 / 6)));
		bary.push_back(Vec3((1.0 / 6), (4.0 / 6), (1.0 / 6)));

		bary.push_back(Vec3((5.0 / 12), (2.0 / 12), (5.0 / 12)));
		bary.push_back(Vec3((2.0 / 12), (5.0 / 12), (5.0 / 12)));
		bary.push_back(Vec3((5.0 / 12), (5.0 / 12), (2.0 / 12)));

		bary.push_back(Vec3((1.0 / 3), (1.0 / 3), (1.0 / 3)));
		break;
	default:
		throw "REALLY ???";
	};

	return bary;
}

Vec3 HW2::projTangent(Vec3 q, Vertex * v){ // project q onto v's tangent plane
	Vec3 tmp =  q - ( v->coords );
	return (v->normal * (tmp*v->normal) )*-1 + q;
}

Vec3 HW2::p(Vec3 b, Vertex** v){
	return  ( v[0]->coords * b.x + v[1]->coords * b.y+  v[2]->coords * b.z); //Vec3(v[0]->coords * b, v[1]->coords * b, v[2]->coords * b);
}

Vec3 HW2::p_star(Vec3 b, Triangle* t){
	Vertex** vs = new Vertex*[3];
	
	vs[ 0 ] = cm->verts[ t->v1i ];
	vs[ 1 ] = cm->verts[ t->v2i ];
	vs[ 2 ] = cm->verts[ t->v3i ];
	
	Vec3 tmp = p(b, vs);
	Vec3 p1 = projTangent( tmp, vs[0] );
	Vec3 p2 = projTangent( tmp, vs[1] );
	Vec3 p3 = projTangent( tmp, vs[2] );

	delete [] vs;
	
	return ( tmp * (1-ALPHA) ) + ( ( p1 * b.x ) + ( p2 * b.y ) + ( p3 * b.z ) ) * ALPHA;
	
	// return Vec3( (b.x * p1.x + b.y * p2.x +  b.z * p3.x),
	// 	     (b.x * p1.y + b.y * p2.y +  b.z * p3.y),
	// 	     (b.x * p1.z + b.y * p2.z +  b.z * p3.z) );
}

Vec3 HW2::findNormal(int v1, int v2, int v3) {
	return (cm->verts[v2]->coords-cm->verts[v1]->coords)^(cm->verts[v3]->coords - cm->verts[v1]->coords);
}

vector<int> HW2::getEdgesM(Mesh* mesh, Triangle* tri){
	vector<int> edgy;

	for (auto ed : mesh->verts[tri->v1i]->edgeList ) 
		if(mesh->edges[ed]->v2i == tri->v2i || mesh->edges[ed]->v1i == tri->v2i )
			edgy.push_back(ed);
	for (auto ed : mesh->verts[tri->v2i]->edgeList ) 
		if(mesh->edges[ed]->v2i == tri->v3i || mesh->edges[ed]->v1i == tri->v3i )
			edgy.push_back(ed);
	for (auto ed : mesh->verts[tri->v3i]->edgeList ) 
		if(mesh->edges[ed]->v2i == tri->v1i || mesh->edges[ed]->v1i == tri->v1i )
			edgy.push_back(ed);
	return edgy;
}

vector<int> HW2::getEdges(Triangle* tri){
	vector<int> edgy;

	for (auto ed : cm->verts[tri->v1i]->edgeList ) 
		if(cm->edges[ed]->v2i == tri->v2i || cm->edges[ed]->v1i == tri->v2i )
			edgy.push_back(ed);
	for (auto ed : cm->verts[tri->v2i]->edgeList ) 
		if(cm->edges[ed]->v2i == tri->v3i || cm->edges[ed]->v1i == tri->v3i )
			edgy.push_back(ed);
	for (auto ed : cm->verts[tri->v3i]->edgeList ) 
		if(cm->edges[ed]->v2i == tri->v1i || cm->edges[ed]->v1i == tri->v1i )
			edgy.push_back(ed);
	return edgy;
}

vector<Vertex*> HW2::triInterior( vector<Vertex *> pset, Vec3 v1, Vec3 v2, Vec3 v3 ){
	Vec3 x;
	Vec3 y;
	switch(parallel(pset)){
	case 0: 
		// do it on xy plane
		x = Vec3(1,0,0);
		y = Vec3(0,1,0);
		break;

	case 1:
		// do it on yz
		x = Vec3(0,1,0);
		y = Vec3(0,0,1);
		break;
	case 2:
		//do it on xz
		x = Vec3(1,0,0);
		y = Vec3(0,0,1);
		break;
	default:
		throw "OOPS";
		break;
	};

	vector<Vec3> proj;
 	float miny = INF, minx = INF;
	for (int i = 0 ; i < pset.size() ; ++i) {
		proj.push_back( Vec3( (pset[i]->coords*x), (pset[i]->coords*y), 0) );
		if (proj[i].y < miny) 
			miny = proj[i].y;
		if (proj[i].x < minx) 
			minx = proj[i].x;
	}
	Vec3 v1p(v1*x, v1*y, 0);
	if (v1p.x < minx) 
		minx = v1p.x;
	if (v1p.y < miny) 
		minx = v1p.y;

	Vec3 v2p(v2*x, v2*y, 0);
	if (v2p.x < minx) 
		minx = v2p.x;
	if (v2p.y < miny) 
		minx = v2p.y;
	
	Vec3 v3p(v3*x, v3*y, 0);	
	if (v3p.x < minx) 
		minx = v3p.x;
	if (v3p.y < miny) 
		minx = v3p.y;
	
	
	Vec3 ori = Vec3(minx, miny, 0);
	for (int i = 0 ; i < pset.size() ; ++i) 
		proj[i] = proj[i] - ori;
	v1p = v1p - ori;
	v2p = v2p - ori;
	v3p = v3p - ori;

	vector<Vertex *> inside;
	for (int i = 0 ; i < pset.size(); ++i)
		if(isInside(proj[i], v1p, v2p, v3p )) 
			inside.push_back(pset[i]);
	return inside;
}

bool HW2::isInside(Vec3 point, Vec3 v1, Vec3 v2, Vec3 v3){ // vec2
	bool inside = 1;
	Vec3 mid = (v1 + v2 + v3)/3;
	
	Vec3 l = (v2-v1); // vec2
	Vec3 normLine(l.x, -(l.x/l.y), 0);
	float mag1 = ((v1 - point) * normLine) / pow(normLine.magnitude(), 2);
	float mag2 = ((v1 - mid  ) * normLine) / pow(normLine.magnitude(), 2);

	if (mag1 == 0) 
		return false;

	
	inside *= ((mag2 > 0) == (mag1 > 0));
	
	l = (v3-v2); // vec2
	normLine = Vec3(l.x, -(l.x/l.y), 0);
	mag1 = ((v2 - point) * normLine) / pow(normLine.magnitude(), 2);
	mag2 = ((v2 - mid) * normLine) / pow(normLine.magnitude(), 2);
	if (mag1 == 0) 
		return false;

	inside *= ((mag2 > 0) == (mag1 > 0));
	
	l = (v1 - v3); // vec2
	normLine = Vec3(l.x, -(l.x/l.y), 0);
	mag1 = ((v2 - point) * normLine) / pow(normLine.magnitude(), 2);
	mag2 = ((v2 - mid) * normLine) / pow(normLine.magnitude(), 2);
	if (mag1 == 0) 
		return false;

	inside *= ((mag2 > 0) == (mag1 > 0));

	return inside;
}

vector<Vec3i> HW2::polygonTriangulation(Mesh* mesh, vector<int> hull){
	vector<Vec3i> triVerts;
	
	float minAngle;
	int idx;
	float aci;
	while (hull.size() > 3) {
		idx = 1;
		minAngle = INF;
		for (int i = 0 ; i < hull.size(); ++i) {
			if (i == 0){
				aci = angleInDeg( mesh->coord( hull[ i+1 ] )            - mesh->coord( hull[ i ] ),
						  mesh->coord( hull[  hull.size()-1 ] ) - mesh->coord( hull[ i ] ) ) ;
			}else if(i == hull.size()-1){
				aci = angleInDeg( mesh->coord( hull[ 0 ] )   - mesh->coord( hull[ i ] ),
						  mesh->coord( hull[ i-1 ] ) - mesh->coord( hull[ i ] ) ) ;

			}else{
				aci = angleInDeg( mesh->coord( hull[ i+1 ] ) - mesh->coord( hull[ i ] ),
						  mesh->coord( hull[ i-1 ] ) - mesh->coord( hull[ i ] ) ) ;
			}

			
			if ( minAngle > aci){
				minAngle = aci;
				idx = i;
			}
		}

		if (idx == 0) {
			triVerts.push_back( Vec3i( hull[0], hull[1], hull[hull.size()-1] ) );
			hull.erase( hull.begin() );
		}else if(idx == hull.size()-1){
			triVerts.push_back( Vec3i( hull[hull.size()-1], hull[0], hull[hull.size()-2]));
			hull.erase( hull.begin() +  hull.size() - 1 );
		}else{
			triVerts.push_back( Vec3i( hull[idx], hull[idx+1], hull[idx-1] ) );
			hull.erase( hull.begin() + idx );
		}
	}
//	triVerts.push_back(Vec3i(hull[0], hull[1], hull[2]));
	return triVerts;
}

float HW2::angleInDeg(Vec3 v1, Vec3 v2){
	return angleInGrad(v1,v2) * 180 / M_PI;
}

double HW2::angleInGrad(Vec3 v1, Vec3 v2){
	return acos( ( v1 * v2 ) / ( v1.magnitude() * v2.magnitude()  ) );
}

void HW2::copyVertex(Mesh *msh, Mesh* sub, int idx){
	msh->verts[idx]->coords = sub->verts[idx]->coords;
	msh->verts[idx]->normal = sub->verts[idx]->normal;
	
	msh->verts[idx]->vertList = sub->verts[idx]->vertList;
	msh->verts[idx]->triList  = sub->verts[idx]->triList;
	msh->verts[idx]->edgeList = sub->verts[idx]->edgeList;
	
	msh->verts[idx]->color = sub->verts[idx]->color;	
};

void HW2::triangulate(Mesh* sub, vector<Vertex*> pset){

	// triangulate the given convex hull for pset
	vector<int> hull = grahamScan(pset);
	vector<Vec3i> nTris = polygonTriangulation(sub, hull);
	
	// add triangles inside each new triangles
	vector<int> edgessIdx;		
	vector<Edge*> edgess;
	vector<Vertex*> inside;
	vector<bool> isHandl;
	for (int i = 0 ; i < 4*sub->verts.size(); ++i) 
		isHandl.push_back(false);
	
	for (int j = 0; j < nTris.size(); ++j) {
		inside = triInterior(pset,
				     sub->verts[nTris[j].x]->coords,
				     sub->verts[nTris[j].y]->coords,
				     sub->verts[nTris[j].z]->coords );
		// CREATE FACES

		vector<int> tmp;		
		tmp = recurTri(sub, inside,
			       sub->verts[nTris[j].x],
			       sub->verts[nTris[j].y],
			       sub->verts[nTris[j].z],
			       isHandl);

		
		for (int i = 0 ; i < tmp.size(); ++i)
			edgessIdx.push_back(tmp[i]);
	}
	
	for (int i = 0 ; i < edgessIdx.size(); ++i)
		edgess.push_back(sub->edges[edgessIdx[i]]);
	
	// DELAUNIZE TRIS
	delanuay(sub, edgess); // INEFFICIENT
}

vector<int> HW2::uniqeSet(vector<int> set ){
	bool flag;
	int curr;
	do {
		for (int i = 0; i < set.size(); ++i) {
			curr = set[i];
			flag = false;
			for (int j = i+1; j < set.size(); ++j) {
				if (curr == set[j]) {
					set.erase(set.begin() + j);
					flag = true;
					break;
				}
			}
			if (flag) {
				break;
			}
		}
	}while(flag);

	return set;
}



vector<int> HW2::recurTri(Mesh* mesh, vector<Vertex*> remaining, Vertex* v1, Vertex* v2, Vertex* v3, vector<bool> &isHandl){

	if (remaining.size() == 0) {
		mesh->addTriangle(v1->idx, v2->idx, v3->idx);
		vector<int> tmp = getEdgesM(mesh, mesh->tris[mesh->tris.size()-1]);
		vector<int> tmpeepeepee;
		if (!isHandl[tmp[0]]) {
			tmpeepeepee.push_back(tmp[0]);
			isHandl[tmp[0]] = true;
		}
		if (!isHandl[tmp[1]]) {
			tmpeepeepee.push_back(tmp[1]);
			isHandl[tmp[0]] = true;
		}
		if (!isHandl[tmp[2]]) {
			tmpeepeepee.push_back(tmp[2]);
			isHandl[tmp[0]] = true;
		}
		
		return tmpeepeepee;
	}

	Vertex* nCorner = remaining[0];
	remaining.erase(remaining.begin()+0);

	vector<Vertex*> r1 = triInterior(remaining, v1->coords, v2->coords, nCorner->coords);
	vector<Vertex*> r2 = triInterior(remaining, v2->coords, v3->coords, nCorner->coords);
	vector<Vertex*> r3 = triInterior(remaining, v3->coords, v1->coords, nCorner->coords);

	vector<int> tmpeepee;
	vector<int> tmpee = recurTri( mesh, r1, v1, v2, nCorner, isHandl);
	for (int i = 0; i < tmpee.size(); ++i)
		if ( !isHandl[tmpee[i]] ) {
			tmpeepee.push_back(tmpee[i]);
			isHandl[tmpee[i]] = true;
		}

	tmpee = recurTri( mesh, r2, v2, v3, nCorner, isHandl);
	for (int i = 0; i < tmpee.size(); ++i) 
		if (!isHandl[tmpee[i]] ) {
			tmpeepee.push_back(tmpee[i]);
			isHandl[tmpee[i]] = true;
		}


	tmpee = recurTri( mesh, r3, v3, v1, nCorner, isHandl);
	for (int i = 0; i < tmpee.size(); ++i) 
		if (!isHandl[tmpee[i]] ) {
			tmpeepee.push_back(tmpee[i]);
			isHandl[tmpee[i]] = true;
		}
	return tmpeepee;
	
	////// AFTER R1 AND R2 INITIALIZATIONS
	///
	// while(flag){
	// 	flag = false;
	// 	for (int i = 0; i < remaining.size(); ++i){
	// 		for (int j = 0; j < r1.size(); ++j){
	// 			if (r1[j] == remaining[i]) {	
	// 				remaining.erase(remaining.begin() + i);
	// 				flag = true;
	// 				break;
	// 			}
	// 		}
	// 		if (flag)
	// 			break;
	// 	}
	// }

	// while(flag){
	// 	flag = false;
	// 	for (int i = 0; i < remaining.size(); ++i){
	// 		for (int j = 0; j < r2.size(); ++j){
	// 			if (r2[j] == remaining[i]) {	
	// 				remaining.erase(remaining.begin() + i);
	// 				flag = true;
	// 				break;
	// 			}
	// 		}
	// 		if (flag)
	// 			break;
	// 	}
	// }	
}

void HW2::delanuay(Mesh* mesh, vector<Edge*> edgess){

	bool isDelanuay;
	do{
		isDelanuay = true;
		for (auto ed : edgess){
			int indicator = ed->isNotSingle();
			switch(indicator){
			case 1:
			case 2:
				continue;
			case 0:
				break;
			default:
				throw "BUT HOW ???";
			}

		
			int v3, v4;
			if      ( mesh->tris[ed->t1i]->v1i != ed->v1i && mesh->tris[ed->t1i]->v1i != ed->v2i ){
				v3 = mesh->tris[ed->t1i]->v1i;
			}else if ( mesh->tris[ed->t1i]->v2i != ed->v1i && mesh->tris[ed->t1i]->v2i != ed->v2i ){
				v3 = mesh->tris[ed->t1i]->v2i;
			}else{//  if ( mesh->tris[ed->t1i]->v3i != ed->v1i && mesh->tris[ed->t1i]->v3i != ed->v2i ){
				v3 = mesh->tris[ed->t1i]->v3i;
			}// }

			if ( mesh->tris[ed->t2i]->v1i != ed->v1i && mesh->tris[ed->t2i]->v1i != ed->v2i ){
				v4 = mesh->tris[ed->t2i]->v1i;
			}else if ( mesh->tris[ed->t2i]->v2i != ed->v1i && mesh->tris[ed->t2i]->v2i != ed->v2i ){
				v4 = mesh->tris[ed->t2i]->v2i;
			}else{//  if ( mesh->tris[ed->t2i]->v3i != ed->v1i && mesh->tris[ed->t2i]->v3i != ed->v2i ){
				v4 = mesh->tris[ed->t2i]->v3i;
			}// }


			Line line1(ed->v1i, v3);

			if (!line1.isOnRight(mesh->verts[v4]->coords) )
				continue;


			//if (180 <= angleInDeg(mesh->coord(v3) - mesh->coord(ed->v1i),
			//		      mesh->coord(v3) - mesh->coord(ed->v1i) )){ 
			//	continue;
			//}
			if ( minAngleForFlip(mesh->coord(ed->v1i), mesh->coord(ed->v2i), mesh->coord(v3), mesh->coord(v4))
			     <
			     minAngleForFlip(mesh->coord(v3), mesh->coord(v4), mesh->coord(ed->v2i), mesh->coord(ed->v1i))){ // FLIPEED
				flip( mesh, mesh->tris[ed->t1i],
				      mesh->tris[ed->t2i],
				      ed, ed->v1i,
				      ed->v2i,
				      v3,
				      v4
				      );
				isDelanuay = false;
			}
		}
	}while(!isDelanuay);

	return;
}

void HW2::flip(Mesh* mesh, Triangle* t1, Triangle* t2, Edge * e, int v1, int v2, int v3, int v4){

	t1->v1i = v3;
	t1->v2i = v2;
	t1->v3i = v4;

	t2->v1i = v4;
	t2->v2i = v1;
	t2->v3i = v3;	

	e->v1i = v3;
	e->v2i = v4;

	// VERT UPDATE

	// CHANGE FROM V1	
	for (int i = 0 ; i < mesh->verts[v1]->vertList.size() ; ++i) 
		if ( mesh->verts[v1]->vertList[i] == v2) {
			mesh->verts[v1]->vertList.erase(mesh->verts[v1]->vertList.begin() + i);
			break;
		}

	for (int i = 0 ; i < mesh->verts[v1]->edgeList.size() ; ++i) 
		if ( mesh->edges[mesh->verts[v1]->edgeList[i]]->idx == e->idx) {
			mesh->verts[v1]->edgeList.erase(mesh->verts[v1]->edgeList.begin() + i);
			break;
		}
	for (int i = 0 ; i < mesh->verts[v1]->triList.size() ; ++i) 
		if ( mesh->tris[mesh->verts[v1]->triList[i]]->idx == t2->idx) {
			mesh->verts[v1]->triList.erase(mesh->verts[v1]->triList.begin() + i);
			break;
		}

	// CHANGE FROM V2
	for (int i = 0 ; i < mesh->verts[v2]->vertList.size() ; ++i) 
		if ( mesh->verts[v2]->vertList[i] == v1) {
			mesh->verts[v2]->vertList.erase(mesh->verts[v2]->vertList.begin() + i);
			break;
		}

	for (int i = 0 ; i < mesh->verts[v2]->edgeList.size() ; ++i) 
		if ( mesh->edges[mesh->verts[v2]->edgeList[i]]->idx == e->idx) {
			mesh->verts[v2]->edgeList.erase(mesh->verts[v2]->edgeList.begin() + i);
			break;
		}
	for (int i = 0 ; i < mesh->verts[v2]->triList.size() ; ++i) 
		if ( mesh->tris[mesh->verts[v2]->triList[i]]->idx == t1->idx) {
			mesh->verts[v2]->triList.erase(mesh->verts[v2]->triList.begin() + i);
			break;
		}

	// V3 & V4 vertLish UPDATE
	mesh->makeVertsNeighbor(v3, v4);

	// CHANGE FROM V3
	mesh->verts[v3]->triList.push_back(t2->idx);
	mesh->verts[v3]->edgeList.push_back(e->idx);

	// CHANGE FROM V4	
	mesh->verts[v4]->triList.push_back(t1->idx);
	mesh->verts[v4]->edgeList.push_back(e->idx);


}

float HW2::minAngleForFlip(Vec3 v1, Vec3 v2, Vec3 v3, Vec3 v4){
	vector<float> arr;

	Vec3 m = v2 - v1;

	Vec3 f1 = v3 - v2;
	Vec3 f2 = v4 - v2;

	Vec3 b1 = v3 - v1;
	Vec3 b2 = v4 - v1;
	
	// front
	arr.push_back( acos( ( (m*-1)  * f1      ) / ( m.magnitude()  * f1.magnitude() ) ) * 180 / M_PI );
	arr.push_back( acos( ( f2      * (m*-1 ) ) / ( m.magnitude()  * f2.magnitude() ) ) * 180 / M_PI );
	// sides
	arr.push_back( acos( ( (f1*-1) * (b1*-1) ) / ( b1.magnitude() * f1.magnitude() ) ) * 180 / M_PI );
	arr.push_back( acos( ( (b2*-1) * (f2*-1) ) / ( b2.magnitude() * f2.magnitude() ) ) * 180 / M_PI );
	// back
	arr.push_back( acos( ( b1      * m       ) / ( b1.magnitude() * m.magnitude()  ) ) * 180 / M_PI );
	arr.push_back( acos( ( m       * b2      ) / ( b2.magnitude() * m.magnitude()  ) ) * 180 / M_PI );


	float min = INF;
	for (int i = 0; i < arr.size(); ++i)
		if (arr[i] < min) 
			min = arr[i];
	
	return min;
}


int HW2::parallel(vector<Vertex*> pset){
	
	// Vec3 before = pset[0]->coords;
	// for (int i = 1; i < pset.size(); ++i) {
	// 	if (pset[i]->coords.x != before.x) 
	// 		return 1; // PARALLEL TO YZ
	// 	if (pset[i]->coords.y != before.y) 
	// 		return 2; // PARALLEL TO XZ

	// 	before = pset[i]->coords;
	// }

	return 0; // PARALLEL TO XY or NO PARALLELITY
}

vector<int> HW2::grahamScan(vector<Vertex *> pset){
	
	// CREATE BARYCENTRIC COORDINATES FOR THE VERTECIES
	vector<Vec3> proj;
	vector<float> angles;
	vector<int> order;


	Vec3 x;
	Vec3 y;

	int tmpeepee = parallel(pset);
	// PROJECT TO 2D
	switch(tmpeepee){
	case 0: 
		// do it on xy plane
		x = Vec3(1,0,0);
		y = Vec3(0,1,0);
		break;

	case 1:
		// do it on yz
		x = Vec3(0,1,0);
		y = Vec3(0,0,1);
		break;
	case 2:
		//do it on xz
		x = Vec3(1,0,0);
		y = Vec3(0,0,1);
		break;
	default:
		throw "OOPS";
		break;
	};
	
	float tmp;
	float miny = INF, minx = INF;
	int curPy = -1  , curPx = -1;
	for (int i = 0 ; i < pset.size() ; ++i) {
		proj.push_back( Vec3((pset[i]->coords*x), (pset[i]->coords*y), 0) );
		if (proj[i].y < miny) {
			miny = proj[i].y;
			curPy = i;
		}
		if (proj[i].x < minx) {
			minx = proj[i].x;
			curPx = i;
		}		
	}
	Vec3 ori = Vec3(minx, miny, 0);
	for (int i = 0 ; i < pset.size() ; ++i) 
		proj[i] = proj[i] - ori;



	for (auto coord : proj) {
		if (coord == proj[curPy]) {
			continue;
		}
		Vec3 curEdge = (coord - proj[curPy]);
		angles.push_back( acos( ( curEdge * (x * -1) ) /
					( curEdge.magnitude() * x.magnitude()) ) * 180 / M_PI);
	}

	vector<float> tmpr;
	while(angles.size() != 0){
		float max = -1;
		int id = -1;
		
		for (int i = 0; i < angles.size(); ++i) {
			if ( angles[i] > max ) {
				max  = angles[i];
				id = i;
			}
		}
		tmpr.push_back(max);
		order.push_back(id);

		angles.erase(angles.begin() + id);
	}

	vector<int> result;
	result.push_back(curPy);

	for (int i = 0 ; i < order.size(); ++i) {
		result.push_back(order[i]);

		if (result.size() == 2) 
			continue;
		bool flag = true;
		float aci;
		while(flag){
			if (result.size() == 2)
			{
				break;
			}
			Line line1(proj[result[result.size() - 1]], proj[result[result.size() - 2]]);
			if(line1.isOnRight(proj[result[result.size() - 1]]))
				result.erase( result.begin() + result.size()-2 );
			else
				flag = false;
		}
	}

	result.erase(result.begin() + result.size()-1);
	result.erase(result.begin() + result.size()-1);
	return result;
}

#endif
 
