#include "Graphics.h"

/*@note: Draws at location (x,y) with specified colour (pColor) at framebuffer (pFrameBuffer)
@param x: X coordinate where to draw the pixel
@param y: Y coordinate where to draw the pixel
@param pColor: Colour to be used
@param pFrameBuffer: Framebuffer where the pixel will be drawn*/
void Draw(int32_t x, int32_t y, GL::SRGB& pColor, float* pFrameBuffer)
{
	assert(pFrameBuffer);
	pFrameBuffer[3 * (y * GL::nWindowWidth + x) + 0] = pColor.r;
	pFrameBuffer[3 * (y * GL::nWindowWidth + x) + 1] = pColor.g;
	pFrameBuffer[3 * (y * GL::nWindowWidth + x) + 2] = pColor.b;
}

/*@note: Sets the mode to be used for rendering, either wireframe or solid filled
@param PolygonMode: The kind of mode needed*/
void SetPolygonMode(EPolygonMode PolygonMode)
{
	GL::PolygonMode = PolygonMode;
}

/*@note: Transforms normalized device coordinates to screen coordinates, does not take depth into account at the moment!
@param x: X coordinate of the viewport rectangle
@param y: Y coordinate of the viewport rectangle
@param nWidth: Viewport width
@param nHeight: Viewport height
@param NDCCoordinates: Normalized device coordinates to be transformed to screen coordinates*/
void ViewportTransform(int32_t x, int32_t y, int32_t nWidth, int32_t nHeight, MVector4& NDCCoordinates, MVector4& Screen)
{
	Screen.x = (nWidth >> 1)		*	NDCCoordinates.x + (x + (nWidth >> 1));
	Screen.y = (nHeight >> 1)	*	NDCCoordinates.y + (y + (nHeight >> 1));
	printf("Screen coords: %f, %f\n", Screen.x, Screen.y);
}