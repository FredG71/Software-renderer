#pragma once
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <gl/freeglut.h>
#include <math.h>
#include <cstdint>
#include <float.h>
#include <cstdio>
#include <cassert>

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
	/* Colour in RGB colour space, where r is Red, g is Green, b is Blue - S is a prefix*/
	struct SRGB{
		float r;
		float g;
		float b;

		SRGB(float __r, float __g, float __b)
			: r(__r), g(__g), b(__b)
		{}
	};
};

void Draw(int32_t x, int32_t y, GL::SRGB& pColor, float* pFrameBuffer);
void Render();