#pragma once
#include <cstdint>
#include "Graphics.h"
#include "Math.h"

/*Primitive types, used in the Draw procedure.*/
enum class EPrimitive{
	// Draws a line
	EPrimitive_Type_Line,
	// Draws a point
	EPrimitive_Type_Point,
	// Draws a triangle
	EPrimitive_Type_Triangle
};

/* Holds an array of vectors */
struct Buffer{
	// Default constructor serves no purpose
	Buffer() = delete;

	/*@param nNum: Number of vectors 3 for a triangle etc.*/
	explicit Buffer(uint32_t nNum)
	:	pBufferAttribute(new MVector4[sizeof(MVector4)* nNum]),
		nNumVectors(nNum)
	{}

	/*
	@param x: X value
	@param y: Y value
	@param z: Z value
	@param w: W value
	@param nIndex: Which vector should we set*/
	inline void SetAttribute(float x, float y, float z, float w, uint32_t nIndex)
	{
		assert(pBufferAttribute && &pBufferAttribute[nIndex]);
		// Ensure we're not writing past the framebuffer!
		assert(	x < GL::nWindowHeight && x < GL::nWindowWidth &&
				y < GL::nWindowHeight && y < GL::nWindowWidth);
		pBufferAttribute[nIndex].x = x;
		pBufferAttribute[nIndex].y = y;
		pBufferAttribute[nIndex].z = z;
		pBufferAttribute[nIndex].w = w;
	}

	~Buffer()
	{
		delete[] pBufferAttribute;
	}

	MVector4* pBufferAttribute;
	uint32_t nNumVectors;
};

void DDALine(int32_t x1, int32_t y1, int32_t x2, int32_t y2, GL::SRGB& Colour);
void DrawPrimitive(Buffer& VertexBuffer, GL::SRGB& Colour, EPrimitive PrimitiveType);
void Rasterize(const MVector4& A, const MVector4& B, const MVector4& C, GL::SRGB& pColour);
void Circle(int32_t x0, int32_t y0, int32_t nRadius);