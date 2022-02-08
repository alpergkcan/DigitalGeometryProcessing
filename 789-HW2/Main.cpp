#define _CRT_SECURE_NO_DEPRECATE
#include <inttypes.h>
#define HAVE_INT8_T
#include <cmath>
#include <xlocnum>
#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>

#include "HW2.h"

int main(int argc, char ** argv)
{
  // Scene Initialization
	HWND window = SoWin::init(argv[0]);
	SoWinExaminerViewer * viewer = new SoWinExaminerViewer(window);
	SoSeparator * root = new SoSeparator;
	root->ref();

	// stuff to be drawn on screen must be added to the root
	Mesh* mesh = new Mesh;
	Painter* painter = new Painter();

	mesh->loadOff((char*) ASSET);
	HW2* sol = new HW2(mesh);

	
#ifdef FINAL
	if (sol->exit)
		goto quit;

	if (sol->paint)
		painter->paint = true;
	else
		painter->paint = false;

	if (sol->astar) {
		root->addChild(painter->getShapeSep(sol->cm));
		root->addChild(painter->getSphereSep(sol->cm, 0, 0, 1));
	}
	else if (sol->subdiv) {
		root->addChild(painter->getShapeSep(sol->nm));
		root->addChild(painter->getSphereSep(sol->nm, 0, 0, 1));
	}
	else {
		root->addChild(painter->getShapeSep(sol->cm));
	}
#endif

#ifndef FINAL
	sol->alpha = ALPHA;
#endif // !FINAL


	// Visualization
#ifdef SUB_DIV_
	root->addChild( painter->getShapeSep(sol->cm) );
#endif
#if SUB_DIV_3 ||  SUB_DIV_4 || SUB_DIV_P 
	root->addChild( painter->getShapeSep(sol->nm) );
#endif


#ifdef VISUALIZE
	root->addChild( painter->getSphereSep(sol->cm, 0, 0, 1) );
#endif
	
	viewer->setSize(SbVec2s(1360, 768));
	viewer->setSceneGraph(root);
	viewer->show();
	SoWin::show(window);
	SoWin::mainLoop();
quit:
	delete viewer;
	root->unref();


	for (int i=0; i < sol->nm->verts.size(); ++i)
		delete sol->nm->verts[i];
	delete sol;
	return 0;
}
