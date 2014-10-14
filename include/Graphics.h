#pragma once
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <gl/freeglut.h>
#include <math.h>
#include <cstdint>
#include <float.h>
#include <cstdio>
#include <cassert>
#include "Math.h"

enum class EPolygonMode{
	// Draws in wireframe mode, i.e. no rasterization occurs.
	EPolygon_Wireframe,
	// Rasterizes the triangle
	EPolygon_Filled
};

namespace GL{
	// Window height - changed on the reshape callback!
	extern int32_t nWindowHeight;
	// Window width - changed on the reshape callback!
	extern int32_t nWindowWidth;
	// Size of framebuffer"
	extern uint32_t nSize;
	// Contains the data to be sent to the framebuffer
	extern float* pFrameBuffer;
	// Old size, used to reallocate the framebuffer as needed
	extern uint32_t nOldSize;
	// Polygon mode - either wireframe or filled
	extern EPolygonMode PolygonMode;
	// Screen coordinates transformed from NDC space - see ::ViewportTransform
	extern MVector4 ScreenCoordinates;
	/* Colour in RGB colour space, where r is Red, g is Green, b is Blue - S is a prefix*/
	struct SRGB{
		float r;
		float g;
		float b;

		SRGB()
			: r(0), g(0), b(0)
		{}
		SRGB(float __r, float __g, float __b)
			: r(__r), g(__g), b(__b)
		{}
		SRGB& operator* (float f)
		{
			r = r * f;
			g = g * f;
			b = b * f;
			return *this;
		}

		SRGB& operator*(SRGB& inColor)
		{
			r = r * inColor.r;
			g = g * inColor.g;
			b = b * inColor.b;
			return *this;
		}
		 SRGB& operator - ( SRGB& inColor)
		{
			r = r - inColor.r;
			g = g - inColor.g;
			b = b - inColor.b;
			return *this;
		}

		 SRGB& operator + (SRGB& inColor)
		 {
			 r = r + inColor.r;
			 g = g + inColor.g;
			 b = b + inColor.b;
			 return *this;
		 }
	};
};

void Draw(int32_t x, int32_t y, GL::SRGB& pColor, float* pFrameBuffer);
void Render();
void SetPolygonMode(EPolygonMode PolygonMode);