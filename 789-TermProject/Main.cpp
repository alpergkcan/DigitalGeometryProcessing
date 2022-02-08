#include "pch.h"

#include "Mesh.h"
#include "input.h"
#include "ARAP.h"
#include "Painter.h"


//#define ASSET  "./mesh/stick_start.obj"
//#define CHANGE_TO "./mesh/stick_end.obj"

#define ASSET  "./mesh/low1.obj"
#define CHANGE_TO "./mesh/low2.obj"


#define WIDTH  1024
#define HEIGHT 640




int main(int argc, char** argv)
{
	// Scene Initialization
	HWND* window = new HWND();
	*window = SoWin::init(argv[0]);

	SoWinExaminerViewer* viewer = new SoWinExaminerViewer(*window);
	root = new SoSeparator;
	root->ref();

	// stuff to be drawn on screen must be added to the root
	start = new Mesh;
	end   = new Mesh;
	painter = new Painter();

	//mesh->loadOff( (string)ASSET );
	//mesh->loadObj( (const char*)ASSET);

	start->loadObj((const char*)ASSET);
	end->loadObj((const char*)CHANGE_TO);

	/*start->push_left(SCREEN_SIZE);
	  end->push_right(SCREEN_SIZE);*/



	for (size_t i = 0; i < 200 ; i++)
	{
		size_t tmp = std::rand() % end->verts.size();
		start->samples.push_back(tmp);
		end->samples.push_back(  tmp);
	}

	//startSep = painter->getShapeSep(start);
	//endSep = painter->getShapeSep(end);

	ARAP* sol = new ARAP(start, end);
	
	int a;
	std::cout << "How many interpolated meshes do you want to create?\n  -> ";
	std::cin >> a;
	std::cout << "\n";


	//for (int i = 0; i < sol->newMeshes.size(); i++) {
	//	for (Vertex* vert : sol->newMeshes[i]->verts)
	//		vert->coord /= 1e2;
	//}
	sol->create_x_interpolations(a);
	sol->place_evenly();
	sol->draw(root, painter);
	startSamp = painter->getSphereSep(start, 0, 0, 1);
	endSamp = painter->getSphereSep(end, 0, 0, 1);


	//root->addChild(startSep);
	root->addChild(startSamp);
	//root->addChild(endSep);
	root->addChild(endSamp);

	//for (size_t i = 0; i < a+2; i++)
	//	root->addChild(painter->getShapeSep(sol->newMeshes[i]));

	//sol->create_interpolation(1.0f);
	//root->addChild(painter->getShapeSep(sol->newMeshes[1]));



	//Mesh* slm = new Mesh;
	//for (size_t i = 0; i < sol->n; i++)
	//{
	//	Matrix3f* M = sol->Ms[sol->cm->verts[i]->triList[0]];
	//	Vector3f* T = sol->Ts[sol->cm->verts[i]->triList[0]];

	//	vec3 coord = (*M) * (sol->cm->verts[i]->coord) + (*T);
	//	slm->addVertex(coord);
	//}
	//for (size_t i = 0; i < sol->size; i++)
	//	slm->addTriangle(sol->cm->tris[i]->v1i, sol->cm->tris[i]->v2i, sol->cm->tris[i]->v3i);
	//root->addChild(painter->getShapeSep(slm));

	//Matrix3f deneme = MatrixXf::Random(3,3);
	//Matrix3f deneme;
	//deneme << 1, 2, 3, 5, 4, 8, 9, 5, 2;
	//std::cout << deneme << "\nNORM: " << deneme.norm() << "\n";
	////Quat q = QH::to_quat_recip(deneme).normalized();
	//Quat qq = QH::to_quat(deneme).normalized();
	////std::cout << q << " " << q.norm() << "\n\n";
	//std::cout << qq << " " << qq.norm() << "\n\n";
	////Matrix3f deneme_c = QH::to_mat(q);
	//Matrix3f deneme_b = QH::to_mat(qq);
	////std::cout << deneme_c << "\n\n";
	//std::cout << deneme_b << "\n\n";



	//JacobiSVD<MatrixXf> svd(deneme, ComputeFullU | ComputeFullV);

	//Mat3 P = svd.matrixU();
	//Mat3 Q = svd.matrixV();
	//std::cout << P << "\n\n" << Q << "\n\n";

	//Vec3 D_vec = svd.singularValues();
	//Mat3 D;
	//D << D_vec[0], 0, 0, 0, D_vec[1], 0, 0, 0, D_vec[2];
	//std::cout << D << "\n\n";
	//sol->create_interpolation(1);
	//root->addChild(painter->getShapeSep(sol->newMeshes[sol->newMeshes.size()-1]));

	//Mat3 I = mat3::Identity();
	//std::cout << I;



	SbBool(*fcnPtr)(void*, MSG*) { &inputCallback};
	viewer->setEventCallback((SoWinRenderAreaEventCB*)(fcnPtr));
	viewer->setSize( SbVec2s(WIDTH, HEIGHT) );
	viewer->setSceneGraph(root);
	viewer->show();
	SoWin::show(*window);
	SoWin::mainLoop();
quit:
	delete viewer;
	root->unref();
	delete sol;	
	return 0;
}
