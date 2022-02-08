#include "pch.h"
#include "ARAP.h"
#include "Painter.h"
#define SCREEN_SIZE 25


void ARAP::place_evenly() {
	int tt = newMeshes.size();
	float curtis_conner = 1; float max_conner = FLT_MIN;

	for (size_t i = 0; i < tt; i++)
	{
		curtis_conner = newMeshes[i]->maxx - newMeshes[i]->minx;
		if (curtis_conner > max_conner)
			max_conner = curtis_conner;
	}

	float screen_size = (tt)* curtis_conner;
	//float screen_size = 25;

	float move = -screen_size - newMeshes[0]->minx;
	float step = (screen_size - newMeshes[newMeshes.size() - 1]->maxx - move) / (tt - 1);
	for (int i = 0; i < tt; i++) {
		for (Vertex* vert : newMeshes[i]->verts) {
			vert->coord.x() += move;
		}
		newMeshes[i]->minx += move;
		newMeshes[i]->maxx += move;
		move += step;
	}
		
}
void ARAP::draw(SoSeparator* root, Painter* paint) {
	for (int i = 0; i < newMeshes.size(); i++)
		root->addChild(paint->getShapeSep(newMeshes[i]));
}

void ARAP::ARAP_main() {
	Ms.reserve(size);
	Ts.reserve(size);
	Rs.reserve(size);
	Ss.reserve(size);

	alpha = sqrt(pow(max(cm->maxy, tm->maxy) - min(cm->miny, tm->miny), 2) +
				 pow(max(cm->maxx, tm->maxx) - min(cm->minx, tm->minx), 2) +
				 pow(max(cm->maxz, tm->maxz) - min(cm->minz, tm->minz), 2) );
	alpha = 1 / alpha;
	alpha = alpha * alpha;

	
	mat3 *M;
	mat3 P, D, Q, *R, *S;	
	Vec3* T;
	vector<mat3*> svd ;
	for (size_t i = 0; i < size; ++i) {
		Vec3 *start_v1, * start_v2, * start_v3;
		Vec3* end_v1, * end_v2, * end_v3;



		Vec3 start_v4, end_v4, cr, com;
		start_v1 = &(cm->verts[cm->tris[i]->v1i]->coord);
		start_v2 = &(cm->verts[cm->tris[i]->v2i]->coord);
		start_v3 = &(cm->verts[cm->tris[i]->v3i]->coord);
		cr = (*start_v2 - *start_v1).cross(*start_v3 - *start_v2);
		//cr   = (*start_v2 - *start_v1).cross(*start_v3 - *start_v1);
		com = (*start_v1 + *start_v2 + *start_v3) / 3;
		start_v4 = com + cr / sqrt(cr.norm());

		end_v1 = &(tm->verts[tm->tris[i]->v1i]->coord);
		end_v2 = &(tm->verts[tm->tris[i]->v2i]->coord);
		end_v3 = &(tm->verts[tm->tris[i]->v3i]->coord);
		cr = (*end_v2 - *end_v1).cross(*end_v3 - *end_v2);
		//cr = (*end_v2 - *end_v1).cross(*end_v3 - *end_v1);
		com = (*end_v1 + *end_v2 + *end_v3) / 3;
		end_v4 = com + 	cr / sqrt(cr.norm());



		vec3 c1 = *end_v1 - end_v4;
		vec3 c2 = *end_v2 - end_v4;
		vec3 c3 = *end_v3 - end_v4;
		//std::cout << c1 << "\n\n" << c2 << "\n\n" << c3 << "\n";
		mat3 V;
		V.col(0) = c1;
		V.col(1) = c2;
		V.col(2) = c3;
		//std::cout << V << "\n\n";

		c1 = *start_v1 - start_v4;
		c2 = *start_v2 - start_v4;
		c3 = *start_v3 - start_v4;
		//std::cout << c1 << "\n\n" << c2 << "\n\n" << c3 << "\n";
		mat3 U;
		U.col(0) = c1;
		U.col(1) = c2;
		U.col(2) = c3;
		//std::cout << U << "\n\n";



		M = new mat3( V * U.inverse() );
		Ms.push_back( M );

		T = new vec3( (end_v4) - (*M) * (start_v4) );
		Ts.push_back( T );


		JacobiSVD<MatrixXf> svd(*M, ComputeFullU | ComputeFullV);

		P = svd.matrixU();
		Q = svd.matrixV().conjugate().transpose();
		Vec3 D_vec = svd.singularValues();
		D << D_vec[0],		  0,		0,
					0, D_vec[1],		0,
					0,		  0, D_vec[2];
		
		R = new mat3( P*Q );                   
		Rs.push_back( R );
		S = new mat3( Q.transpose()*D*Q );
		Ss.push_back( S );
	}
}

	
void ARAP::create_interpolation(float t) {
	Mesh* newMesh = new Mesh();
	interp_Ms.reserve(size);
	interp_Ts.reserve(size);
	Eigen::Quaternion<float> zero_rot;

	zero_rot = {1,0,0,0};
	float* errors = new float[size];

	mat3 Mt, Rt, St;
	vec3 Tt;
	mat3 I = mat3::Identity();

	for (size_t i = 0; i < size; ++i) {
		Quaternion<float> rot_quat, interp_quat;

		rot_quat    = QH::to_quat( *(Rs[i]) );

		//interp_quat = zero_rot.slerp( t , rot_quat);
		interp_quat = zero_rot.slerp(t, rot_quat);
		// interp results
		Rt = QH::to_mat(interp_quat);
		St = (I * (1 - t) + *Ss[i] * t);
		Mt = Rt * St;
		Tt = *(Ts[i]) * t;
        interp_Ms.push_back( new Mat3( Mt ) );
		interp_Ts.push_back( new Vec3( Tt ) );

		
		// ERROR CALCULATION
		float A = cm->tris[i]->getArea(cm);
		float M_related = ((*Ms[i]) - *interp_Ms[i]).norm(); // CHECK IF CONJUGATE !!!
		float T_related = ((*Ts[i]) - *interp_Ts[i]).norm();
		float error = A * (M_related * M_related + alpha * T_related * T_related);

		errors[i] = error;
	}

	for (Vertex* vert : cm->verts){
		float min_error = FLT_MAX;
		size_t min_id = -1;
		for (size_t id : vert->triList) {
			if ( errors[id] < min_error) {
				min_id = id;
				min_error = errors[id];
			}
		}

		vec3 coor = interpolate_coord(vert->coord, min_id);
		//vec3 coor = interpolate_coord(vert->coord, vert->triList[0]);

		newMesh->addVertex(coor);
	}
	for (Triangle* face: cm->tris)
		newMesh->addTriangle(face->v1i, face->v2i, face->v3i);

	for (size_t i = 0; i < size; i++)
	{
		delete interp_Ms[i];
		delete interp_Ts[i];
	}
	delete[]errors;
	interp_Ms.clear();
	interp_Ts.clear();

	newMeshes.push_back(newMesh);	
}

void ARAP::create_x_interpolations(int x) {
	float step = 1.0f/(x+1.0f);
	float time = step;
	newMeshes.push_back(cm);

	for (int i = 0; i < x; time += step, i++) 
		create_interpolation(time);
	newMeshes.push_back(tm);
}
void ARAP::normal_fix() {
	for (size_t i = 0; i < newMeshes.size(); i++)
	{
		for (Vertex* vert : newMeshes[i]->verts)
		{
			if (vert->coord.norm() < FLT_EPSILON)
			{
				
			}
		}
	}
}

//vector<Tetrahedron*> ARAP::tetrahedralize(Mesh* msh){
//	vector<Tetrahedron*> array;
//	for (Triangle* face: msh->tris)
//		array.push_back(new Tetrahedron(msh, face));
//	return array;
//}



vec3 ARAP::interpolate_coord(const vec3 &coord, size_t transform_idx) {
	Matrix3f* Mt = interp_Ms[transform_idx];
	Vector3f* Tt = interp_Ts[transform_idx];

	return  *Mt * coord + *Tt;
}


/// UTILITY
Vec3 ARAP::findNormal(int v1, int v2, int v3) {
	return (cm->verts[v2]->coord - cm->verts[v1]->coord).cross(cm->verts[v3]->coord - cm->verts[v1]->coord ).normalized();
}

vector<int> ARAP::getEdgesM(Mesh* mesh, Triangle* tri){
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

vector<int> ARAP::getEdges(Triangle* tri){
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

float ARAP::angleInDeg(Vec3 v1, Vec3 v2){
	return angleInGrad(v1,v2) * 180 / M_PI;
}

double ARAP::angleInGrad(Vec3 v1, Vec3 v2){
	return acos( v1.dot(v2) / (v1.norm() * v2.norm()) ) ;
}

void ARAP::copyVertex(Mesh *msh, Mesh* sub, int idx){
	msh->verts[idx]->coord = sub->verts[idx]->coord;
	msh->verts[idx]->normal = sub->verts[idx]->normal;
	
	msh->verts[idx]->vertList = sub->verts[idx]->vertList;
	msh->verts[idx]->triList  = sub->verts[idx]->triList;
	msh->verts[idx]->edgeList = sub->verts[idx]->edgeList;
	
	msh->verts[idx]->color = sub->verts[idx]->color;	
};

