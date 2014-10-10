#pragma once
#include "CPU.h"
#include <math.h>
#include <smmintrin.h>
#include <cassert>


/* 4 dimensional vector class, currently used as MVector2 and MVector3 have not been implemented yet
@todo: Implement cross product, normalization, operator overloading (vector addition, scalar multiplication etc)*/
class _MM_ALIGN16 MVector4{
public:
	inline MVector4() :
		nSSE4(GetSSE4Support())
	{
		if (nSSE4)
			mVector = _mm_setzero_ps();
		else
			x = y = z = 0;
	}
	inline MVector4(float x, float y, float z, float w) :
		nSSE4(GetSSE4Support())
	{
		if (nSSE4)
			mVector = _mm_set_ps(w, z, y, x);
		else
		{
			this->x = x;
			this->y = y;
			this->z = z;
			this->w = w;
		}
	}
	explicit MVector4(float a) :
		nSSE4(GetSSE4Support())
	{
		if (nSSE4)
			mVector = _mm_set_ps(a, a, a, a);
		else
			x = y = z = w = a;
	}
	inline void Dot4(MVector4& Other)
	{
		assert(nSSE4);
		mVector = _mm_dp_ps(mVector, Other.mVector, 0xFF);
	}
	inline void Dot4(MVector4& A, MVector4& Other)
	{
		assert(nSSE4);
		mVector = _mm_dp_ps(A.mVector, Other.mVector, 0xFF);
	}
	inline float fDot4(MVector4& A)
	{
		return (x * A.x + y * A.y + z * A.z + w * A.w);
	}
	inline float fDot4(MVector4& A, MVector4& Other)
	{
		return (A.x* Other.x + A.y * Other.y + A.z * Other.z + A.w * Other.w);
	}

	/*Copy constructor/assignment and move constructor/assignment not implemented yet!*/
	MVector4(const MVector4&) = delete;
	MVector4(MVector4&&) = delete;
	MVector4& operator=(const MVector4&) = delete;
	MVector4& operator=(MVector4&&) = delete;

	const int32_t nSSE4;
	union{
		struct{
			float x;
			float y;
			float z;
			float w;
		};
		__m128 mVector;
	};
};