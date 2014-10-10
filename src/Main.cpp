#include "Primitive.h"
#include <time.h>

/*
*	Author: Frederic Garnier.
*	Language: C++
*	Created: 9/10/2014
*	Modified: 10/10/2014
*
*
*	Currently supports drawing lines, among other things see debugCheckerBoard(), no support for linear interlation, matrices, 2 dimensional vectors
*	3 dimensional vectors etc, this should work on other platforms, however it hasn't been tested
*	
*	Libraries: OpenGL, Freeglut.
*	Compiler: Microsoft Visual C++ Compiler Nov 2013 CTP
*
*	Platform tested: Windows 7 x64
*
*
*	BUGS:
*	A bug can occur while drawing a line, such that the line will not go from the start of the monitor to the end. More than one triangle and more than one
*	line is not supported at the moment, and can cause undefined results.
*/
int32_t GL::nWindowHeight = 256;
int32_t GL::nWindowWidth = 256;
uint32_t GL::nSize = GL::nWindowWidth * GL::nWindowHeight;
float* GL::pFrameBuffer;
uint32_t GL::nOldSize = GL::nSize;
EPolygonMode GL::PolygonMode = EPolygonMode::EPolygon_Filled;


Buffer VertexBuffer = Buffer(3);

/*@note: Draws a procedural checkerboard for debugging purposes*/
void debugCheckerBoard()
{
	GL::SRGB FColour = GL::SRGB(0.f, 0.f, 0.f);
	for (uint32_t x = 0; x < GL::nWindowWidth; ++x)
	{
		for (uint32_t y = 0; y < GL::nWindowHeight; ++y)
		{
			GLubyte byCheck = (GLubyte)(((y & 0x10) == 0) ^ ((x & 0x10) == 0));
			FColour.g = byCheck;
			FColour.b = byCheck;
			Draw(x, y, FColour, GL::pFrameBuffer);
		}
	}
}

/*@note: Main loop, this is where any drawing is to be done
@todo: Allocate pFrameBuffer initially once, and only reallocate when window is resized, allocating and deallocating like so, ought to be very very slow.*/
void Render()
{
	// Adjust array size as window is resized
	GL::nSize = GL::nWindowWidth * GL::nWindowHeight;
	GL::pFrameBuffer = new float[GL::nSize * 3];
	DrawPrimitive(VertexBuffer, GL::SRGB(1, 1, 0), EPrimitive::EPrimitive_Type_Triangle);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDrawPixels(GL::nWindowWidth, GL::nWindowHeight, GL_RGB, GL_FLOAT, GL::pFrameBuffer);
	glutSwapBuffers();
	delete[] GL::pFrameBuffer;
}


/*@note: Re size callback, called whenever the window needs to be resized, the window's width and height 
(stored in the GL namespace) are updated accordingly
@param nWidth: Width of the window
@param nHeight: Height of the window*/
void OnReshape(int32_t nWidth, int32_t nHeight)
{
	GL::nWindowHeight = nHeight;
	GL::nWindowWidth = nWidth;
}

void InitBuffer()
{
	VertexBuffer.SetAttribute(20, 20, 0, 0, 0);
	VertexBuffer.SetAttribute(GL::nWindowHeight - 1, GL::nWindowHeight - 1, 0, 0, 1);
	VertexBuffer.SetAttribute(GL::nWindowHeight - 1, 20, 0, 0, 2);
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);
	glutInitWindowSize(GL::nWindowWidth, GL::nWindowHeight);
	glutCreateWindow("Frederic Garnier - Software rendering");
	InitBuffer();
	glutDisplayFunc(Render);
	glutReshapeFunc(OnReshape);
	glEnable(GL_DEPTH_TEST);
	glClearColor(0.f, 0.f, 0.f, 1.f);
	glutMainLoop();
}