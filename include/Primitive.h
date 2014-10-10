#pragma once
#include <cstdint>
#include <math.h>
#include "Graphics.h"
#include "Vector.h"

/*Primitive types, used in the Draw procedure.*/
enum class EPrimitive{
	// Draws a line
	EPrimitive_Type_Line,
	// Draws a point
	EPrimitive_Type_Point,
	// Draws a triangle
	EPrimitive_Type_Triangle
};

void DDALine(int32_t x1, int32_t y1, int32_t x2, int32_t y2, GL::SRGB& Colour);
void DrawPrimitive(MVector4& Position, MVector4& End, GL::SRGB& Colour, EPrimitive PrimitiveType);