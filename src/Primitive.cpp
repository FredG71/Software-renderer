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

/*@note: Rasterizes a triangle with input vectors A, B, C and colour pColour
@param A: Coordinate of one of the point of the triangle
@param B: Coordinate of one of the point of the triangle
@param C: Coordinate of one of the point of the triangle
@param pColour: Colour of the resulting triangle*/
void Rasterize(const MVector4& A, const MVector4& B, const MVector4& C, GL::SRGB& pColour)
{
	float ya = A.y;
	float yb = B.y;
	float yc = C.y;

	float xa = A.x;
	float xb = B.x;
	float xc = C.x;


	register int32_t nMinx = (int32_t)min(min(xa, xb), xc);
	register int32_t nMaxx = (int32_t)max(max(xa, xb), xc);
	register int32_t nMiny = (int32_t)min(min(ya, ya), yc);
	register int32_t nMaxy = (int32_t)max(max(ya, yb), yc);

	for (int32_t y = nMiny; y < nMaxy; y++)
	{
		for (int32_t x = nMinx; x < nMaxx; x++)
		{
			if (((xa - xb) * (y - ya) - (ya - yb) * (x - xa) > 0) &&
				((xb - xc) * (y - yb) - (yb - yc) * (x - xb) > 0) &&
				((xc - xa) * (y - yc) - (yc - ya) * (x - xc) > 0)
				)
			{
				GL::pFrameBuffer[3 * (y * GL::nWindowWidth + x) + 0] = pColour.r;
				GL::pFrameBuffer[3 * (y * GL::nWindowWidth + x) + 1] = pColour.g;
				GL::pFrameBuffer[3 * (y * GL::nWindowWidth + x) + 2] = pColour.b;
			}
		}
	}
}

/*@note: Draws a circle using the circle midpoint algorith, see wikipedia page for implementation
@param x0: X position of the individual points
@param y0: Y position of the individual points
@param nRadius: Radius of the circle*/
void Circle(int32_t x0, int32_t y0, int32_t nRadius)
{
	int32_t x = nRadius;
	int32_t y = 0;
	int nRadiusError = 1 - x;

	while (x >= y)
	{
		Draw(x + x0, y + y0, GL::SRGB(1, 0, 0), GL::pFrameBuffer);
		Draw(y + x0, x + y0, GL::SRGB(1, 0, 0), GL::pFrameBuffer);
		Draw(-x + x0, y + y0, GL::SRGB(1, 0, 0), GL::pFrameBuffer);
		Draw(-y + x0, x + y0, GL::SRGB(1, 0, 0), GL::pFrameBuffer);
		Draw(-x + x0, -y + y0, GL::SRGB(1, 0, 0), GL::pFrameBuffer);
		Draw(-y + x0, -x + y0, GL::SRGB(1, 0, 0), GL::pFrameBuffer);
		Draw(x + x0, -y + y0, GL::SRGB(1, 0, 0), GL::pFrameBuffer);
		Draw(y + x0, -x + y0, GL::SRGB(1, 0, 0), GL::pFrameBuffer);
		y++;

		if (nRadiusError < 0)
			nRadiusError += 2 * y + 1;
		else
		{
			x--;
			nRadiusError += 2 * (y - x + 1);
		}
	}
}

/*@note: Draws a primitive. Lines and triangles and points are supported, however more than one line and more than one triangle hasn't been tested yet,
multiple points works fine, no colour interpolation occurs.
@param VertexBuffer: A buffer that holds the location of where to render (As 4 dimensional vector at the moment).
@param Colour: Colour of the line
@param PrimitiveType: Type of primitive*/
void DrawPrimitive(Buffer& VertexBuffer, GL::SRGB& Colour, EPrimitive PrimitiveType)
{
	switch (PrimitiveType)
	{
	case EPrimitive::EPrimitive_Type_Line:
		assert(VertexBuffer.nNumVectors % 2 == 0);
		DDALine(VertexBuffer.pBufferAttribute[0].x, VertexBuffer.pBufferAttribute[0].y,
			VertexBuffer.pBufferAttribute[1].x, VertexBuffer.pBufferAttribute[1].y, Colour);
		break;
	case EPrimitive::EPrimitive_Type_Point:
		if (VertexBuffer.nNumVectors == 1)
		{
			Draw(VertexBuffer.pBufferAttribute[0].x, VertexBuffer.pBufferAttribute[0].y, Colour, GL::pFrameBuffer);
			break;
		}
		for (uint32_t nIndex = 0; nIndex < VertexBuffer.nNumVectors; nIndex++)
			Draw(VertexBuffer.pBufferAttribute[nIndex].x, VertexBuffer.pBufferAttribute[nIndex].y, Colour, GL::pFrameBuffer);
		break;
	case EPrimitive::EPrimitive_Type_Triangle:
		assert(VertexBuffer.nNumVectors % 3 == 0);
		switch (GL::PolygonMode)
		{
		case EPolygonMode::EPolygon_Filled:
			Rasterize(VertexBuffer.pBufferAttribute[0], VertexBuffer.pBufferAttribute[1], VertexBuffer.pBufferAttribute[2], Colour);
			break;
		case EPolygonMode::EPolygon_Wireframe:
			DDALine(VertexBuffer.pBufferAttribute[0].x, VertexBuffer.pBufferAttribute[0].y,
				VertexBuffer.pBufferAttribute[1].x, VertexBuffer.pBufferAttribute[1].y, Colour);
			DDALine(VertexBuffer.pBufferAttribute[1].x, VertexBuffer.pBufferAttribute[1].y,
				VertexBuffer.pBufferAttribute[2].x, VertexBuffer.pBufferAttribute[2].y, Colour);
			DDALine(VertexBuffer.pBufferAttribute[2].x, VertexBuffer.pBufferAttribute[2].y,
				VertexBuffer.pBufferAttribute[0].x, VertexBuffer.pBufferAttribute[0].y, Colour);
			break;
		}
		break;
	}
}