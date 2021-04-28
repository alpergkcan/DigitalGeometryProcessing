#define _CRT_SECURE_NO_DEPRECATE
#include <inttypes.h>
#define HAVE_INT8_T
#include <cmath>
#include <xlocnum>

#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>

#include "HW1.h"

int main(int, char ** argv)
{
  // Scene Initialization
	HWND window = SoWin::init(argv[0]);
	SoWinExaminerViewer * viewer = new SoWinExaminerViewer(window);
	SoSeparator * root = new SoSeparator;
	root->ref();

	// stuff to be drawn on screen must be added to the root
	Mesh* mesh = new Mesh();
	Painter* painter = new Painter();
	
	// HW1 Main
	mesh->loadOff((char*)"assets/woman.off");
	HW1* sol = new HW1(mesh);

	// Visualization
	root->addChild( painter->getShapeSep(mesh) );
	root->addChild( painter->getSphereSep(mesh, 0, 0, 1) );
	viewer->setSize(SbVec2s(640, 480));
	viewer->setSceneGraph(root);
	viewer->show();
	SoWin::show(window);
	SoWin::mainLoop();
	delete viewer;
	root->unref();
	return 0;
}
