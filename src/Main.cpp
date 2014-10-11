#include "Primitive.h"
#include <time.h>

/*
*	Author: Frederic Garnier.
*	Language: C++
*	Created: 9/10/2014
*	Modified: 11/10/2014
*
*
*	Currently supports drawing lines, among other things see debugCheckerBoard(), no support for linear interlation, 2 dimensional vectors
*	3 dimensional vectors etc, this should work on other platforms, however it hasn't been tested
*	
*	Libraries: OpenGL, Freeglut.
*	Compiler: Microsoft Visual C++ Compiler Nov 2013 CTP
*
*	Platform tested: Windows 7 x64
*
*
*	BUGS:
*	Rotation around the Z axis can cause issues, which will be resolved when going to 3D, stretching the window can cause a crash, due to going out of the
*	framebuffer's bounds etc.s
*/
int32_t GL::nWindowHeight = 512;
int32_t GL::nWindowWidth = 512;
uint32_t GL::nSize = GL::nWindowWidth * GL::nWindowHeight;
float* GL::pFrameBuffer;
uint32_t GL::nOldSize = GL::nSize;
EPolygonMode GL::PolygonMode = EPolygonMode::EPolygon_Filled;
MVector4 GL::ScreenCoordinates = MVector4(0);

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
	VertexBuffer.SetAttribute(-1, -1, 0, 0, 0);
	VertexBuffer.SetAttribute(0, 1, 0, 0, 1);
	VertexBuffer.SetAttribute(1, -1, 0, 0, 2);
	MMatrix4x4 pScaleMatrix = ScaleMatrix(pScaleMatrix, 0.5, 0.5, 0.5);
	MMatrix4x4 pRotationMatrix = RotationMatrixY(pRotationMatrix, 1.1);
	MMatrix4x4 MatrixProduct = pScaleMatrix * pRotationMatrix;
	for (int nIndex = 0; nIndex < VertexBuffer.nNumVectors; nIndex++)
	{
		VertexBuffer.pBufferAttribute[nIndex] = VertexBuffer.pBufferAttribute[nIndex] * MatrixProduct;
	}
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