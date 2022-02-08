#pragma once
#include "pch.h"

#include "Mesh.h"
#include "Painter.h"
//#include "ARAP.h"


// mouse states
#define MOUSE_ROAM 0 
#define MOUSE_STILL 1 

#define MOUSE_BUTTON1_CLICK 2
#define MOUSE_BUTTON1_PUSH 3 
#define MOUSE_BUTTON1_ROAM 4
#define MOUSE_BUTTON1_RELEASE 5 

#define MOUSE_BUTTON2 6 

#define MOUSE_BUTTON3_CLICK 7
#define MOUSE_BUTTON3_PUSH 8
#define MOUSE_BUTTON3_ROAM 9
#define MOUSE_BUTTON3_RELEASE 10

#define MOUSE_RELEASE 11
// keyboard states
#define KEYBOARD_FIRST 12
#define KEYBOARD_STICKY 13
#define KEYBOARD_STICKY_RELEASE 14
#define KEYBOARD_RELEASE 15
#define KEYBOARD_CLICK 16

#define KEYBOARD_WAITING 17


// messages
#define RAW_KEYBOARD_FIRST 256
#define RAW_KEYBOARD_PRESS 258
#define RAW_KEYBOARD_RELEASE 257

#define RAW_MOUSE_WEIRD 132
#define RAW_MOUSE_ROAM 512
#define RAW_MOUSE_BUTTON1_PRESS 513
#define RAW_MOUSE_BUTTON1_RELEASE 514
#define RAW_MOUSE_BUTTON2 516
#define RAW_MOUSE_BUTTON3_PRESS 519
#define RAW_MOUSE_BUTTON3_RELEASE 520
#define RAW_MOUSE_RELEASE 533

Mesh* mesh;
SoSeparator* root;
SoSeparator* startSep;
SoSeparator* endSep;

SoSeparator* startSamp;
SoSeparator* endSamp;

Mesh* start, * end;
Painter* painter;

	//void init(Mesh* msh, ARAP* arap, SoSeparator* sep1, SoSeparator* sep2, Painter* paint);

	SbBool changed = FALSE;

	uint state = MOUSE_ROAM;
	uint prev_state = MOUSE_ROAM;
	char kb_input = 64;
	SbBool kb_updated = FALSE;

	//void changeMesh(string name) {
	//	mesh->changeObj(name.c_str());
	//	root->removeChild(mshSep);
	//	mshSep = painter->getShapeSep(mesh);
	//	root->addChild(mshSep);
	//}

	void printMSG(MSG* message) {
		std::cout << "Message: " << message->message << std::endl
			<< "wParam:  " << message->wParam << std::endl
			<< "lParam:  " << message->lParam << std::endl
			<< "time:    " << message->time << std::endl
			<< "coord:   " << message->pt.x << " " << message->pt.y << std::endl << std::endl;
	}
	 uint translateMSG(MSG* message, uint last_state) {
		switch (message->message)
		{
		case RAW_KEYBOARD_FIRST:
			if ((last_state == KEYBOARD_CLICK) | (last_state == KEYBOARD_STICKY))
				return KEYBOARD_STICKY;
			else
				return KEYBOARD_FIRST;
			break;
		case RAW_KEYBOARD_PRESS:
			if (last_state == KEYBOARD_STICKY)
			{
				return KEYBOARD_STICKY;
			}
			else if (last_state == KEYBOARD_FIRST) {
				return KEYBOARD_CLICK;
			}
			else
				std::cout << "errorPress\n";
			break;
		case RAW_KEYBOARD_RELEASE:
			if (last_state == KEYBOARD_CLICK)
			{
				return KEYBOARD_RELEASE;
			}
			else if (last_state == KEYBOARD_STICKY) {
				return KEYBOARD_STICKY_RELEASE;
			}
			else
			{
				std::cout << "errorRelease\n";
			}
			break;



		case RAW_MOUSE_BUTTON1_PRESS:
			//std::cout << "mb1 pressed" << std::endl;
			return MOUSE_BUTTON1_PUSH;
			break;
		case RAW_MOUSE_BUTTON1_RELEASE:
			if (last_state == MOUSE_BUTTON1_PUSH)
				return MOUSE_BUTTON1_CLICK;
			else if (last_state == MOUSE_BUTTON1_ROAM)
				return MOUSE_BUTTON1_RELEASE;
			else
				std::cout << "error1\n" << std::endl;
			break;
		case RAW_MOUSE_ROAM:
			if ((last_state == MOUSE_BUTTON1_PUSH) | (last_state == MOUSE_BUTTON1_ROAM))
			{
				return MOUSE_BUTTON1_ROAM;
			}
			else if ((last_state == MOUSE_BUTTON3_PUSH) | (last_state == MOUSE_BUTTON3_ROAM))
			{
				return MOUSE_BUTTON3_ROAM;
			}
			else {
				return MOUSE_ROAM;
			}
			break;
		case RAW_MOUSE_BUTTON2:
			return MOUSE_BUTTON2;
			break;
		case RAW_MOUSE_BUTTON3_PRESS:
			return MOUSE_BUTTON3_PUSH;
			break;
		case RAW_MOUSE_BUTTON3_RELEASE:
			if (last_state == MOUSE_BUTTON3_PUSH)
				return MOUSE_BUTTON1_CLICK;
			else if (last_state == MOUSE_BUTTON3_ROAM)
				return MOUSE_BUTTON3_RELEASE;
			else
				std::cout << "error3\n" << std::endl;
			break;
		default:
			return last_state;
			break;
		}
	}

	 SbBool inputCallback(void* userData, MSG* anyevent)
	{
		SbBool handled = FALSE;

		state = translateMSG(anyevent, prev_state);

		//SoWinRenderArea* myRenderArea = (SoWinRenderArea*)userData;
		//SoMouseButtonEvent* myMouseButtonEvent;

		SbVec3f vec;

		//// print
		//printMSG(anyevent);
		//std::cout << "state: " << state << std::endl;
		////

		switch (state) {

		case KEYBOARD_CLICK: // mb1
			kb_input = anyevent->wParam;
			kb_updated = true;
			handled = TRUE;
			break;


		default:
			// handled = FALSE;
			break;
		}

		if (kb_updated) {
			switch (kb_input)
			{
			case 65:
			case 97:
				break;
			case 1111111111:
				if (changed)
				{
					changed = false;
					//changeMesh((string)ASSET);
					root->removeChild(endSep);
					root->removeChild(endSamp);

					startSep = painter->getShapeSep(start);
					startSamp = painter->getSphereSep(start, 0, 0, 1);
					root->addChild(startSep);
					root->addChild(startSamp);

				}
				else
				{
					changed = true;
					//changeMesh((string)CHANGE_TO);
					root->removeChild(startSep);
					root->removeChild(startSamp);

					endSep  = painter->getShapeSep(end);
					endSamp = painter->getSphereSep(end, 0, 0, 1);
					root->addChild(endSep);
					root->addChild(endSamp);

				}
				break;
			default:
				break;
			}
			kb_updated = false;
		}


		prev_state = state;
		return handled;
	}