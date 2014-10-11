#pragma once
#include "CPU.h"
#include <math.h>
#include <smmintrin.h>
#include <cassert>

struct MMatrix4x4;
class _MM_ALIGN16 MVector4;


__m128 SSE_TransformVector(MVector4& Vector, MMatrix4x4& Matrix);
#define SSE_SetVector( A, Index ) _mm_shuffle_ps( A, A, (( Index ) | ((  Index ) << 2 ) | ( Index ) << 4 |  ( Index ) << 6 ))
inline __m128 SSE_MulVector(__m128 Vector, __m128 Matrix);
inline __m128 SSE_AddVector(__m128 VectorA, __m128 VectorB);

/* 4 dimensional vector class, currently used as MVector2 and MVector3 have not been implemented yet
@todo: Implement cross product, normalization, operator overloading (vector addition, scalar multiplication etc)*/
class _MM_ALIGN16 MVector4{
public:
	MVector4& operator= (MVector4& Other)
	{
		this->x = Other.x;
		this->y = Other.y;
		this->z = Other.z;
		this->w = Other.w;
		return *this;
	}
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
	inline MVector4& operator* (MVector4& Other)
	{
		x *= Other.x;
		y *= Other.y;
		z *= Other.z;
		w *= Other.w;
		return *this;
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

	inline MVector4& operator *(MMatrix4x4& Other)
	{
		static MVector4 Temp = MVector4(0);
		Temp.mVector = SSE_TransformVector(*this, Other);
		return Temp;
	}
	const int32_t nSSE4;
	union{
		struct{
			float x;
			float y;
			float z;
			float w;
		};
		float pVectorArray[4];
		__m128 mVector;
	};
};


/* Matrix class, no SSE support at this moment!*/
struct MMatrix4x4{
public:
	MMatrix4x4()
		:	m00(1), m01(0), m02(0), m03(0),
			m10(0), m11(1), m12(0), m13(0),
			m20(0), m21(0), m22(1), m23(0),
			m30(0), m31(0), m32(0), m33(1)
	{}
	explicit MMatrix4x4(float f)
		:	m00(f), m01(f), m02(f), m03(f),
			m10(f), m11(f), m12(f), m13(f),
			m20(f), m21(f), m22(f), m23(f),
			m30(f), m31(f), m32(f), m33(f)
	{}

	MMatrix4x4(float __m00, float __m01, float __m02, float __m03,
		float __m10, float __m11, float __m12, float __m13,
		float __m20, float __m21, float __m22, float __m23,
		float __m30, float __m31, float __m32, float __m33
		)
		:	m00(__m00), m01(__m01), m02(__m02), m03(__m03),
			m10(__m10), m11(__m11), m12(__m12), m13(__m13),
			m20(__m20), m21(__m21), m22(__m22), m23(__m23),
			m30(__m30), m31(__m31), m32(__m32), m33(__m33)
	{}
	/*@todo: Figure out how to not return a local variable like so, currently fine since it's static */
	inline MMatrix4x4& operator* (const MMatrix4x4& A)
	{
		static MMatrix4x4 Product = MMatrix4x4(0);
		for (register uint32_t i = 0; i < 4; i++)
		{
			for (register uint32_t j = 0; j < 4; j++)
			{
				for (register uint32_t k = 0; k < 4; k++)
				{
					Product.pMatrixArray[i][j] += pMatrixArray[i][k] * A.pMatrixArray[k][j];
				}
			}
		}
		return Product;
	}

	inline MMatrix4x4& operator* (float fScalar)
	{
		for (register uint32_t i = 0; i < 4; i++)
		{
			for (register uint32_t j = 0; j < 4; j++)
			{
				pMatrixArray[i][j] *= fScalar;
			}
		}
		return *this;
	}

	union{
		struct{
			float	m00, m01, m02, m03,
					m10, m11, m12, m13,
					m20, m21, m22, m23,
					m30, m31, m32, m33;
		};
		float pMatrixArray[4][4];
		__m128 mMatrix[4];
	};
};

/*@param Vector: Vector to be multiplied by the matrix
@param Matrix: Matrix to multiply the vector with
@return: Multiplied vector*/
inline __m128 SSE_MulVector(__m128 Vector, __m128 Matrix)
{
	return _mm_mul_ps(Vector, Matrix);
}

/*@param VectorA: Augend
@param VectorB: Addend
@return: Sum of the two vectors*/
inline __m128 SSE_AddVector(__m128 VectorA, __m128 VectorB)
{
	return _mm_add_ps(VectorA, VectorB);
}

/*@note: Transforms Vector by Matrix
@param Vector: Vector to be transformed
@param Matrix: Target matrix
@return: Transformed vector*/
inline __m128 SSE_TransformVector(MVector4& Vector, MMatrix4x4& Matrix)
{
	__m128 A, B, C, D;
	A = B = C = D = {};

	A = SSE_SetVector(Vector.mVector, 0);
	B = SSE_SetVector(Vector.mVector, 1);
	C = SSE_SetVector(Vector.mVector, 2);
	D = SSE_SetVector(Vector.mVector, 3);

	A = SSE_MulVector(A, Matrix.mMatrix[0]);
	B = SSE_MulVector(B, Matrix.mMatrix[1]);
	C = SSE_MulVector(C, Matrix.mMatrix[2]);
	D = SSE_MulVector(D, Matrix.mMatrix[3]);

	A = SSE_AddVector(A, B);
	C = SSE_AddVector(C, D);
	A = SSE_AddVector(A, C);

	return A;
}

/*@note: Creates a scaling matrix
@param inTarget: the scaling matrix
@param fScaleX: Scale in the x direction
@param fScaleY: Scale in the y direction
@param fScaleZ: Scale in the z direction
@return: Scaling matrix*/
inline MMatrix4x4 ScaleMatrix(MMatrix4x4& inTarget, float fScaleX, float fScaleY, float fScaleZ)
{
	inTarget = MMatrix4x4(	fScaleX, 0, 0, 0,
							0, fScaleY, 0, 0,
							0, 0, fScaleZ, 0,
							0, 0, 0, 1 );
	return inTarget;
}

/*@note: Creates a translation matrix
@param inTarget: the translation matrix
@param fOffsetX: Translation in the x direction
@param fOffsetY: Translation in the y direction
@param fOffsetZ: Translation in the z direction
@return: Translation matrix*/
inline MMatrix4x4 TranslationMatrix(MMatrix4x4& inTarget, float fOffsetX, float fOffsetY, float fOffsetZ)
{
	inTarget = MMatrix4x4(	1, 0, 0, fOffsetX,
							0, 1, 0, fOffsetY,
							0, 0, 1, fOffsetZ,
							0, 0, 0, 1);
	return inTarget;
}

/*@note: Creates a rotation matrix, (Axis is the X axis).
@param inTarget: the rotation matrix
@param fAngle: Angle in RADIANS!
@return: Rotation matrix*/
inline MMatrix4x4 RotationMatrixX(MMatrix4x4& inTarget, float fAngle)
{
	inTarget = MMatrix4x4(	1, 0, 0, 0,
							0, cosf(fAngle), -sinf(fAngle), 0,
							0, sinf(fAngle), cos(fAngle), 0,
							0, 0, 0, 0);
	return inTarget;
}

/*@note: Creates a rotation matrix, (Axis is the Y axis).
@param inTarget: the rotation matrix
@param fAngle: Angle in RADIANS!
@return: Rotation matrix*/
inline MMatrix4x4 RotationMatrixY(MMatrix4x4& inTarget, float fAngle)
{
	inTarget = MMatrix4x4(	cosf(fAngle), 0, sinf(fAngle), 0,
							0, 1, 0, 0,
							-sinf(fAngle), 0, cosf(fAngle), 0,
							0, 0, 0, 0);
	return inTarget;
}

/*@note: Creates a rotation matrix, (Axis is the Z axis).
@param inTarget: the rotation matrix
@param fAngle: Angle in RADIANS!
@return: Rotation matrix*/
inline MMatrix4x4 RotationMatrixZ(MMatrix4x4& inTarget, float fAngle)
{
	inTarget = MMatrix4x4(	cosf(fAngle), -sinf(fAngle), 0, 0,
							sinf(fAngle), cosf(fAngle), 0, 0,
							0, 0, 1, 0,
							0, 0, 0, 0);

	return inTarget;
}

/*@note: Creates a rotation matrix, (Axis is the unit vector Axis).
@param inTarget: the rotation matrix
@param fAngle: Angle in RADIANS!
@param Axis: The axis to rotate about, Axis ought to be a unit vector.
@return: Rotation matrix*/
inline MMatrix4x4 RotationAxisAngle(MMatrix4x4& inTarget, float fAngle, MVector4& Axis)
{
	inTarget = MMatrix4x4(	cosf(fAngle) + powf(Axis.x, 2) * (1 - cosf(fAngle)),
							Axis.x * Axis.y * (1 - cosf(fAngle)) - Axis.z * sinf(fAngle),
							Axis.x * Axis.z * (1 - cosf(fAngle)) + Axis.y * sinf( fAngle),
							Axis.y * Axis.x * (1 - cosf(fAngle)) + Axis.z * sinf(fAngle), 0,
							cosf(fAngle) + pow(Axis.y, 2) * (1 - cosf(fAngle)),
							Axis.y * Axis.z * (1 - cosf(fAngle)) - Axis.x * sinf(fAngle),
							Axis.z * Axis.x * (1 - cosf(fAngle)) - Axis.y * sinf(fAngle), 0,
							Axis.z * Axis.y * (1 - cosf(fAngle)) + Axis.x * sinf(fAngle),
							cosf(fAngle) + powf(Axis.z, 2) * (1 - cosf(fAngle)), 0,
							0, 0, 0, 0);
}