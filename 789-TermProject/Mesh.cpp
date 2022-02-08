#include "pch.h"
#include "Mesh.h"
void Mesh::changeOff(string name) {
	emptyMesh();
	loadOff(name);
}
void Mesh::changeObj(const char* name) {
	emptyMesh();
	loadObj(name);
}

void Mesh::loadObj(const char* name)
{
	std::ifstream file(name);
	std::string line;

	float x, y, z;
	int n = 0;
	int id1, id2, id3;

	minx = FLT_MAX; miny = FLT_MAX; maxx = FLT_MIN; maxy = FLT_MIN;
	maxz = FLT_MIN; minz = FLT_MAX;
	while (std::getline(file, line))
	{
		if (line[0] == '#' || line[0] == 'o' || line[0] == 's')
		{
			continue;
		}
		if (line[0] == 'v')
		{
			std::stringstream ss(line.substr(2));
			ss >> x >> y >> z;
			addVertex(x, y, z);			
		}
		if (line[0] == 'f')
		{
			std::stringstream ss(line.substr(2));
			ss >> id1 >> id2 >> id3;
			addTriangle(id1 - 1, id2 - 1, id3 - 1);
		}
	}
	file.close();
	//CoMToOrigin();
}


void Mesh::loadOff(string name)
{

	std::cout << " LOAD " << std::endl;
	std::ifstream fin;
	fin.open(name);

	char str[334];
	fin >> str;

	int nVerts, nTris, n, i = 0;
	float x, y, z;

	minx = FLT_MAX; miny = FLT_MAX; maxx = FLT_MIN; maxy = FLT_MIN;


	Vector3f sum(0, 0, 0);
	fin >> nVerts >> nTris >> n;
	while (i++ < nVerts)
	{
		fin >> x >> y >> z;
		addVertex(x, y, z);

		//		std::cout << x << " " << y << " " << z << std::endl;
		/*sum = sum + verts[(uint)i - 1]->coord;*/
	}

	int tmp;
	int _x, _y, _z;
	while (fin)
	{
		fin >> tmp;
		fin >> _x >> _y >> _z;

		addTriangle(_x, _y, _z);
	}
	fin.close();
	//CoMToOrigin();

}


void Mesh::addTriangle(int v1, int v2, int v3)
{
	int idx = tris.size();
	tris.push_back(new Triangle(idx, v1, v2, v3));

	//set up structure

	verts[v1]->triList.push_back(idx);
	verts[v2]->triList.push_back(idx);
	verts[v3]->triList.push_back(idx);




	if (!makeVertsNeighbor(v1, v2))
		addEdge(v1, v2);
	else
		edgeUpdate(v1, v2, idx);

	if (!makeVertsNeighbor(v1, v3))
		addEdge(v1, v3);
	else
		edgeUpdate(v1, v3, idx);

	if (!makeVertsNeighbor(v2, v3))
		addEdge(v2, v3);
	else
		edgeUpdate(v2, v3, idx);
}

bool Mesh::makeVertsNeighbor(int v1i, int v2i)
{
	//returns true if v1i already neighbor w/ v2i; false o/w

	for (uint i = 0; i < verts[v1i]->vertList.size(); i++)
		if (verts[v1i]->vertList[i] == v2i)
			return true;

	verts[v1i]->vertList.push_back(v2i);
	verts[v2i]->vertList.push_back(v1i);
	return false;
}
void Mesh::edgeUpdate(int v1, int v2, int idx) {
	int edgeId;
	for (auto ed : verts[v1]->edgeList) {
		if (edges[ed]->v2i == v2 || edges[ed]->v1i == v2) {
			edgeId = ed;
			break;
		}
	}

	int sing = edges[edgeId]->isNotSingle();
	switch (sing) {
	case 0:
		//sedges.push_back(edges[edgeId]);
		//std::cout << "0" << std::endl;
		break;
	case 1:
		edges[edgeId]->t2i = idx;
		break;
	case 2:
		edges[edgeId]->t1i = idx;
		break;
	default:
		throw "REACH TO NONSINGULAR EDGE";
		break;
	};

}

void Mesh::addVertex(float x, float y, float z)
{
	verts.push_back(new Vertex(verts.size(), x, y, z));



	if (x < minx) minx = x;
	else if (x > maxx) maxx = x;

	if (y < miny) miny = y;
	else if (y > maxy) maxy = y;

	if (z < minz) minz = z;
	else if (z > maxz) maxz = z;
}

void Mesh::addVertex(Vec3 nP)
{
	verts.push_back(new Vertex(verts.size(), nP(0), nP(1), nP(2) ));
	if (nP.x() < minx) minx = nP.x();
	else if (nP.x() > maxx) maxx = nP.x();

	if (nP.y() < miny) miny = nP.y();
	else if (nP.y() > maxy) maxy = nP.y();

	if (nP.z() < minz) minz = nP.z();
	else if (nP.z() > maxz) maxz = nP.z();

}

void Mesh::addEdge(int v1, int v2)
{
	int idx = edges.size();
	int t1 = -1, t2 = -1;

	for (auto tri1 : verts[v1]->triList)
		for (auto tri2 : verts[v2]->triList)
			if (tri1 == tri2) {
				t1 = tri1;
				break;
			}

	if (t1 != -1) {
		if (tris[t1]->v1i == v1 && tris[t1]->v3i == v2) {
			t2 = -1;
		}
		else {
			t2 = t1;
			t1 = -1;
		}
	}

	else
		throw "UNREASONABLE AMOUNT OF TRIANGLES AROUND THE EDGE";

	edges.push_back(new Edge(idx, v1, v2, t1, t2));

	verts[v1]->edgeList.push_back(idx);
	verts[v2]->edgeList.push_back(idx);
}

float Mesh::getLength(int v1, int v2) {
	return (verts[v2]->coord - verts[v1]->coord).norm(); // euclidian distance
}

float Mesh::getArea(Triangle* trip) {

	return ((verts[trip->v2i]->coord - verts[trip->v1i]->coord).cross(
		    (verts[trip->v3i]->coord - verts[trip->v1i]->coord)       )).norm() / 2; // cross/2
}

float Edge::getLength(Mesh* msh) { 
	return (msh->verts[v2i]->coord - msh->verts[v1i]->coord).norm(); // euclidian distance
}

float Triangle::getArea(Mesh* msh) {
	return ((msh->verts[v2i]->coord - msh->verts[v1i]->coord).cross(
			(msh->verts[v3i]->coord - msh->verts[v1i]->coord)       ).norm()) / 2; // cross/2
}




int Edge::isNotSingle() {
	return (t1i == -1) * 2 + (t2i == -1) * 1;
}


void  Mesh::push_left(float x) {
	float move_left = -x - minx;
	for (size_t i = 0; i < verts.size(); i++)
		verts[i]->coord.x() += move_left;	
}

void Mesh:: push_right(float x) {
	float move_left = x- maxx;
	for (size_t i = 0; i < verts.size(); i++)
		verts[i]->coord.x() += move_left;
}



void Mesh::emptyMesh() {
	for (luint i = 0; i < verts.size(); ++i) {
		delete verts[i];
	}
	for (luint i = 0; i < edges.size(); ++i) {
		delete edges[i];
	}
	for (luint i = 0; i < tris.size(); ++i) {
		delete tris[i];
	}
	verts.clear();
	edges.clear();
	tris.clear();
	sedges.clear();
	samples.clear();
}



// CONSTRUCTORS
Mesh::Mesh() {
	minx = 0;
	miny = 0;
	minz = 0;
	maxx = 0;
	maxy = 0;
	maxz = 0;
}
Mesh::Mesh(string name) {
	loadOff(name);
}
Mesh::~Mesh() {
	emptyMesh();
}
Vertex::Vertex(int i, float x = 0.0f, float y = 0.0f, float z = 0.0f) {
	idx = i;
	coord = Vec3(x, y, z);
}
Vertex::Vertex(int i, Vertex* vert) {
	idx = i;
	coord = Vec3(vert->coord);
}

Edge::Edge(int id, int v1, int v2) :
	idx(id),
	v1i(v1),
	v2i(v2),
	t1i(-1),
	t2i(-1) {}
Edge::Edge(int id, int v1, int v2, int t1, int t2) :
	idx(id),
	v1i(v1),
	v2i(v2),
	t1i(t1),
	t2i(t2) {}

Triangle::Triangle(int id, int v1, int v2, int v3) :
	idx(id),
	v1i(v1),
	v2i(v2),
	v3i(v3) {};
////


void Mesh::CoMToOrigin()
{
	vec3 sum = { 0, 0, 0 };
	for (Vertex* v : verts)
		sum += v->coord;
	
	sum /= (double)verts.size();

	for (Vertex* v : verts)
		v->coord -= sum;
}