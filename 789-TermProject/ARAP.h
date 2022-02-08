#pragma once
#include "pch.h"
#include "Mesh.h"
#include "Quaternion.h"

class ARAP {
public:
	Mesh* cm; 
	Mesh* tm;
	vector<Mesh*> newMeshes;

	size_t n;
	size_t size;

	float alpha;
	
	vector<mat3*> Ms;
	vector<vec3*> Ts;
	vector<mat3*> Rs;
	vector<mat3*> Ss;

	vector<mat3*> interp_Ms;
	vector<vec3*> interp_Ts;

	

	ARAP(Mesh* start, Mesh* end) : cm(start), tm(end), n(start->verts.size()), size(cm->tris.size()) {
		ARAP_main();
	};

	~ARAP() {
		for (size_t i = 0; i < size; i++)
		{
			delete Ms[i];
			delete Ts[i];
			delete Rs[i];
			delete Ss[i];
		}
	}

 //^^-._________
	void ARAP_main();
	void create_interpolation(float t);
	void create_x_interpolations(int x);
	void normal_fix();
	//vector<Tetrahedron*> tetrahedralize(Mesh* msh);
	vec3 interpolate_coord(const vec3& coord, size_t transform_idx);

	void place_evenly();
	void draw(SoSeparator* root, Painter* paint);

	// utility
	Vec3 findNormal(int v1, int v2, int v3);
	vector<int> getEdges(Triangle* tri);
	vector<int> getEdgesM(Mesh* mesh, Triangle* tri);
	void copyVertex(Mesh* msh, Mesh* sub, int idx);

	float angleInDeg(Vec3 v1, Vec3 v2);
	double angleInGrad(Vec3 v1, Vec3 v2);

	static float* toFloat3(Vec3 r) {
		float* array = new float[3];

		array[0] = r.x();
		array[1] = r.y();
		array[2] = r.z();

		return array;
	}
};

//,,
