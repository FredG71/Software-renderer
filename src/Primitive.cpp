#include "Primitive.h"

/*@note: Draws a line using the differential digital analyzer line algorithm
@param x1: Start x position of the line
@param y1: Start y position of the line
@param x2: End x position of the line
@param y2: End y position of the line*/
void DDALine(
	int32_t x1,
	int32_t y1,
	int32_t x2,
	int32_t y2,
	GL::SRGB& Colour
	)
{
	assert(GL::pFrameBuffer);
	int32_t dx = x2 - x1;
	int32_t dy = y2 - y1;
	int32_t nSteps, nIndex;
	float xIncrement, yIncrement, x = x1, y = y1;

	if (abs(dx) > abs(dy))
		nSteps = abs(dx);
	else
		nSteps = abs(dy);

	xIncrement = dx / (float)nSteps;
	yIncrement = dy / (float)nSteps;

	Draw(roundf(x), roundf(y), Colour, GL::pFrameBuffer);

	for (nIndex = 0; nIndex < nSteps; nIndex++)
	{
		x += xIncrement;
		y += yIncrement;
		Draw(roundf(x), roundf(y), Colour, GL::pFrameBuffer);
	}
}


/*@note: Draws a primitive, line or triangle or point depending upon PrimitiveType, currently uses the DDA algorithm, however
bresenham and possibly Xiaolin Wu's line algorithm ought to be used instead, point and triangle hasn't been implemented yet, and this should be modified
to accomodate for points and triangles (see parameters)
no colour interpolation occurs.
@param Start: Start position (much more suited to lines at the moment)
@param End: End position, same Start
@param Colour: Colour of the line
@param PrimitiveType: Type of primitive*/
void DrawPrimitive(MVector4& Start, MVector4& End, GL::SRGB& Colour, EPrimitive PrimitiveType)
{
	switch (PrimitiveType)
	{
	case EPrimitive::EPrimitive_Type_Line:
		DDALine(Start.x, Start.y, End.x, End.y, Colour);
		break;
	case EPrimitive::EPrimitive_Type_Point:
		break;
	case EPrimitive::EPrimitive_Type_Triangle:
			break;
	default:
		break;
	}
}