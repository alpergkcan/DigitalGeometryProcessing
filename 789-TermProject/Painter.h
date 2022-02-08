#pragma once
#include "pch.h"
#include "Mesh.h"
#include "ARAP.h"



class Painter {
	
public:
	bool paint;
	SoSeparator* getShapeSep(Mesh* mesh);
	SoSeparator* getSphereSep(Mesh* mesh, float deltaX, float deltaY, float scale);
};
