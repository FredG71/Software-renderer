#include "Primitive.h"
#include <time.h>

/*
*	Author: Frederic Garnier.
*	Language: C++
*	Created: 9/10/2014
*	Modified: 14/10/2014
*
*	
*	Libraries: OpenGL, Freeglut.
*	Compiler: Microsoft Visual C++ Compiler Nov 2013 CTP
*
*	Platform tested: Windows 7 x64
*
*
*	BUGS:
*	Multiplying matrix by identify matrix can return a zero matrix, concatenating transforms can cause this issue too. Stretching screen causes issues but that should be
*	solved in the future, matrix functions may have computations errors, need to test those, however the output seems correct.
*/
int32_t GL::nWindowHeight = 768;
int32_t GL::nWindowWidth = 1024;
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
	DrawPrimitive(VertexBuffer, GL::SRGB(1, 0, 1), EPrimitive::EPrimitive_Type_Triangle);
	//debugCheckerBoard();
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
	VertexBuffer.SetAttribute( 0.f, 0.25f, 0.f, 1.f, 0 );
	VertexBuffer.SetAttribute( 0.25f, -0.25f, 0.f, 1.f, 1);
	VertexBuffer.SetAttribute( -0.25, -0.25, 0.f, 1.f, 2 );

	MMatrix4x4 Projection = Perspective(45.f, (float)GL::nWindowWidth / GL::nWindowHeight, 1.f, 100.f);
	MMatrix4x4 View = LookAt(MVector4(2.f, 2.f, -1.f, 0.f), MVector4(0.f, 0.f, 0.f, 0.f), MVector4( 0.f, 1.f, 0.f, 0.f ));

	To3D(0, 0, GL::nWindowWidth, GL::nWindowHeight, View, Projection, VertexBuffer.pBufferAttribute);
	MMatrix4x4 Rotation = RotationMatrixZ(Rotation, 0.3);


	for ( int i = 0; i < 3; i++ )
		VertexBuffer.pBufferAttribute[i] = Rotation * VertexBuffer.pBufferAttribute[i];
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