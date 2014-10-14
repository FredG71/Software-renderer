#pragma once
#define _USE_MATH_DEFINES
#include "CPU.h"
#include <math.h>
#include <smmintrin.h>
#include <cassert>
#include "Debug.h"
struct MMatrix4x4;
class _MM_ALIGN16 MVector4;

void ViewportTransform(int32_t x, int32_t y, int32_t nWidth, int32_t nHeight, MVector4& NDCCoordinates, MVector4& Screen);

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
	inline MVector4(__m128 __mVector)
		: nSSE4(GetSSE4Support()), mVector(__mVector)
	{}
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
	inline float Length()
	{
		return sqrt(powf(x, 2) + powf(y, 2) + powf(z, 2) + powf(w, 2));
	}
	MVector4& operator *(float fIn)
	{
		x *= fIn;
		y *= fIn;
		z *= fIn;
		w *= fIn;
		return *this;
	}

	inline MVector4& Cross(MVector4& A, MVector4& B)
	{
		float __x__ = A.y * B.z - B.y * A.z;
		float __y__ = A.z * B.x - B.z * A.x;
		float __z__ = A.x * B.y - B.x * A.y;

		x = __x__;
		y = __y__;
		z = __z__;

		return *this;
	}

	inline MVector4& Normalize(MVector4& InVector)
	{
		float fLength = InVector.Length();
		mVector = _mm_set_ps(InVector.w, InVector.z, InVector.y, InVector.x);
		x /= fLength;
		y /= fLength;
		z /= fLength;
		return *this;
	}
	inline MVector4& Normalize()
	{
		float fLength = Length();
		x /= fLength;
		y /= fLength;
		z /= fLength;
		return *this;
	}
	MVector4& operator +(MVector4& Other)
	{
		static MVector4 Temp = MVector4(0);
		Temp.x = x + Other.x;
		Temp.y = y + Other.y;
		Temp.z = z + Other.z;
		Temp.w = w + Other.w;
		return Temp;
	}
	MVector4& operator += (MVector4& Other)
	{
		x += Other.x;
		y += Other.y;
		z += Other.z;
		w += Other.w;
		return *this;
	}
	MVector4& operator -(MVector4& Other)
	{
		static MVector4 Temp = MVector4(0);
		Temp.x = x - Other.x;
		Temp.y = y - Other.y;
		Temp.z = z - Other.z;
		Temp.w = w - Other.w;
		return Temp;
	}
	MVector4& operator -= (MVector4& Other)
	{
		x -= Other.x;
		y -= Other.y;
		z -= Other.z;
		w -= Other.w;
		return *this;
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

	MMatrix4x4(	float __m00, float __m01, float __m02, float __m03,
				float __m10, float __m11, float __m12, float __m13,
				float __m20, float __m21, float __m22, float __m23,
				float __m30, float __m31, float __m32, float __m33
		)
		:	m00(__m00), m01(__m01), m02(__m02), m03(__m03),
			m10(__m10), m11(__m11), m12(__m12), m13(__m13),
			m20(__m20), m21(__m21), m22(__m22), m23(__m23),
			m30(__m30), m31(__m31), m32(__m32), m33(__m33)
	{}

	MMatrix4x4& operator* (float fScalar)
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
	return inTarget;
}

/*@note: Normalizes inVector
@param InVector: Vector to be normalized
@return: Normalized vector*/
inline MVector4& Normalize(MVector4& InVector)
{
	static MVector4 Normalized(0);
	float fLength = InVector.Length();
	Normalized.mVector = _mm_set_ps(InVector.w, InVector.z, InVector.y, InVector.x);
	Normalized.x /= fLength;
	Normalized.y /= fLength;
	Normalized.z /= fLength;
	return Normalized;
}

/*@note: Performs a cross product - order A Cross B
@param A: Left hand vector
@param B: Right hand vector
@return: Resulting cross product*/
inline MVector4& Cross(MVector4& A, MVector4& B)
{
	static MVector4 Product(0);
	float __x__ = A.y * B.z - B.y * A.z;
	float __y__ = A.z * B.x - B.z * A.x;
	float __z__ = A.x * B.y - B.x * A.y;

	Product.x = __x__;
	Product.y = __y__;
	Product.z = __z__;

	return Product;
}

/*@note: Returns a look at matrix - RH
@param Position: Camera position
@param Target: Target to look at
@param Up: Up vector
@return: The lookat matrix*/
inline MMatrix4x4 LookAt(MVector4& position, MVector4& target, MVector4& Up)
{
	MVector4 f = Normalize(position - target);
	MVector4 s = Normalize(Cross(f, Up));
	MVector4 u = Cross(s, f);

	MMatrix4x4 ViewMatrix = MMatrix4x4();
	ViewMatrix.pMatrixArray[0][0] = s.x;
	ViewMatrix.pMatrixArray[1][0] = s.y;
	ViewMatrix.pMatrixArray[2][0] = s.z;
	ViewMatrix.pMatrixArray[0][1] = u.x;
	ViewMatrix.pMatrixArray[1][1] = u.y;
	ViewMatrix.pMatrixArray[2][1] = u.z;
	ViewMatrix.pMatrixArray[0][2] = -f.x;
	ViewMatrix.pMatrixArray[1][2] = -f.y;
	ViewMatrix.pMatrixArray[2][2] = -f.z;
	ViewMatrix.pMatrixArray[3][0] = -s.fDot4(s, position);
	ViewMatrix.pMatrixArray[3][1] = -u.fDot4(u, position);
	ViewMatrix.pMatrixArray[3][2] = f.fDot4(f, position);

	return ViewMatrix;

}
/*@note: Returns a perspective matrix - RH
@param FOV: Field of view in degrees
@param aspectRatio: The aspect ratio
@param zNear: z near plane
@param zFar: z far plane
@return: The perspective matrix*/
inline MMatrix4x4& Perspective(float FOV, float aspectRatio, float zNear, float zFar)
{

	float DegToRad = M_PI / 180.f;
	float FovRad = FOV * DegToRad;
	float tanHalfFov = tan(FovRad / 2);
	static MMatrix4x4 Frustrum(1 / (aspectRatio * tanHalfFov), 0, 0, 0,
		0, 1 / tanHalfFov, 0, 0,
		0, 0, -(zFar + zNear) / (zFar - zNear), -1,
		0, 0, -(2 * zFar * zNear) / (zFar - zNear), 0);

	return Frustrum;
}

/*@note: Transforms vector A by Matrix Other (CM)
@param Other: Target matrix
@param A: Target vector
@return: Transformed vector*/
inline MVector4& operator *(MMatrix4x4& Other, MVector4& A)
{
	static MVector4 Product(0);

	Product.x = A.x * Other.pMatrixArray[0][0] +
		A.y * Other.pMatrixArray[0][1] +
		A.z * Other.pMatrixArray[0][2] +
		A.w * Other.pMatrixArray[0][3];

	Product.y = A.x * Other.pMatrixArray[1][0] +
		A.y * Other.pMatrixArray[1][1] +
		A.z * Other.pMatrixArray[1][2] +
		A.w * Other.pMatrixArray[1][3];

	Product.z = A.x * Other.pMatrixArray[2][0] +
		A.y * Other.pMatrixArray[2][1] +
		A.z * Other.pMatrixArray[2][2] +
		A.w * Other.pMatrixArray[2][3];

	Product.w = A.x * Other.pMatrixArray[3][0] +
		A.y * Other.pMatrixArray[3][1] +
		A.z * Other.pMatrixArray[3][2] +
		A.w * Other.pMatrixArray[3][3];
	return Product;
}

/*@note: Multiplies matrix A by matrix B
@param B (left)
@param A: (right)
@return: multiplied matrix
@fix: A bug can occur when multipliying by the identify matrix, or concatenating multiple matrices
i.e. Proj * View * World etc*/
inline MMatrix4x4& operator* (MMatrix4x4& A, MMatrix4x4& B)
{
	static MMatrix4x4 Product = MMatrix4x4(0);

	// COLUMN 0
	Product.pMatrixArray[0][0] = B.pMatrixArray[0][0] * A.pMatrixArray[0][0] +
		B.pMatrixArray[0][1] * A.pMatrixArray[1][0] +
		B.pMatrixArray[0][2] * A.pMatrixArray[2][0] +
		B.pMatrixArray[0][3] * A.pMatrixArray[3][0];

	Product.pMatrixArray[0][1] = B.pMatrixArray[0][0] * A.pMatrixArray[0][1] +
		B.pMatrixArray[0][1] * A.pMatrixArray[1][1] +
		B.pMatrixArray[0][2] * A.pMatrixArray[2][1] +
		B.pMatrixArray[0][3] * A.pMatrixArray[3][1];

	Product.pMatrixArray[0][2] = B.pMatrixArray[0][0] * A.pMatrixArray[0][2] +
		B.pMatrixArray[0][1] * A.pMatrixArray[1][2] +
		B.pMatrixArray[0][2] * A.pMatrixArray[2][2] +
		B.pMatrixArray[0][3] * A.pMatrixArray[3][2];

	Product.pMatrixArray[0][3] = B.pMatrixArray[0][0] * A.pMatrixArray[0][3] +
		B.pMatrixArray[0][1] * A.pMatrixArray[1][3] +
		B.pMatrixArray[0][2] * A.pMatrixArray[2][3] +
		B.pMatrixArray[0][3] * A.pMatrixArray[3][3];


	// COLUMN 1
	Product.pMatrixArray[1][0] = B.pMatrixArray[1][0] * A.pMatrixArray[0][0] +
		B.pMatrixArray[1][1] * A.pMatrixArray[1][0] +
		B.pMatrixArray[1][2] * A.pMatrixArray[2][0] +
		B.pMatrixArray[1][3] * A.pMatrixArray[3][0];

	Product.pMatrixArray[1][1] = B.pMatrixArray[1][0] * A.pMatrixArray[0][1] +
		B.pMatrixArray[1][1] * A.pMatrixArray[1][1] +
		B.pMatrixArray[1][2] * A.pMatrixArray[2][1] +
		B.pMatrixArray[1][3] * A.pMatrixArray[3][1];

	Product.pMatrixArray[1][2] = B.pMatrixArray[1][0] * A.pMatrixArray[0][2] +
		B.pMatrixArray[1][1] * A.pMatrixArray[1][2] +
		B.pMatrixArray[1][2] * A.pMatrixArray[2][2] +
		B.pMatrixArray[1][3] * A.pMatrixArray[3][2];

	Product.pMatrixArray[1][3] = B.pMatrixArray[1][0] * A.pMatrixArray[0][3] +
		B.pMatrixArray[1][1] * A.pMatrixArray[1][3] +
		B.pMatrixArray[1][2] * A.pMatrixArray[2][3] +
		B.pMatrixArray[1][3] * A.pMatrixArray[3][3];

	//COLUMN 2
	Product.pMatrixArray[2][0] = B.pMatrixArray[2][0] * A.pMatrixArray[0][0] +
		B.pMatrixArray[2][1] * A.pMatrixArray[1][0] +
		B.pMatrixArray[2][2] * A.pMatrixArray[2][0] +
		B.pMatrixArray[2][3] * A.pMatrixArray[3][0];

	Product.pMatrixArray[2][1] = B.pMatrixArray[2][0] * A.pMatrixArray[0][1] +
		B.pMatrixArray[2][1] * A.pMatrixArray[1][1] +
		B.pMatrixArray[2][2] * A.pMatrixArray[2][1] +
		B.pMatrixArray[2][3] * A.pMatrixArray[3][1];

	Product.pMatrixArray[2][2] = B.pMatrixArray[2][0] * A.pMatrixArray[0][2] +
		B.pMatrixArray[2][1] * A.pMatrixArray[1][2] +
		B.pMatrixArray[2][2] * A.pMatrixArray[2][2] +
		B.pMatrixArray[2][3] * A.pMatrixArray[3][2];

	Product.pMatrixArray[2][3] = B.pMatrixArray[2][0] * A.pMatrixArray[0][3] +
		B.pMatrixArray[2][1] * A.pMatrixArray[1][3] +
		B.pMatrixArray[2][2] * A.pMatrixArray[2][3] +
		B.pMatrixArray[2][3] * A.pMatrixArray[3][3];

	// COLUMN 3

	Product.pMatrixArray[3][0] = B.pMatrixArray[3][0] * A.pMatrixArray[0][0] +
		B.pMatrixArray[3][1] * A.pMatrixArray[1][0] +
		B.pMatrixArray[3][2] * A.pMatrixArray[2][0] +
		B.pMatrixArray[3][3] * A.pMatrixArray[3][0];

	Product.pMatrixArray[3][1] = B.pMatrixArray[3][0] * A.pMatrixArray[0][1] +
		B.pMatrixArray[3][1] * A.pMatrixArray[1][1] +
		B.pMatrixArray[3][2] * A.pMatrixArray[2][1] +
		B.pMatrixArray[3][3] * A.pMatrixArray[3][1];

	Product.pMatrixArray[3][2] = B.pMatrixArray[3][0] * A.pMatrixArray[0][2] +
		B.pMatrixArray[3][1] * A.pMatrixArray[1][2] +
		B.pMatrixArray[3][2] * A.pMatrixArray[2][2] +
		B.pMatrixArray[3][3] * A.pMatrixArray[3][2];

	Product.pMatrixArray[3][3] = B.pMatrixArray[3][0] * A.pMatrixArray[0][3] +
		B.pMatrixArray[3][1] * A.pMatrixArray[1][3] +
		B.pMatrixArray[3][2] * A.pMatrixArray[2][3] +
		B.pMatrixArray[3][3] * A.pMatrixArray[3][3];

	return Product;
}

/*@note: Transforms vector A by Matrix Other (RM)
@param A: Target vector
@param Other: Target matrix
@return: Transformed vector*/
inline MVector4& operator *(MVector4& A, MMatrix4x4& Other)
{
	static MVector4 Product(0);
	Product.x = A.x * Other.pMatrixArray[0][0] +
		A.y * Other.pMatrixArray[1][0] +
		A.z * Other.pMatrixArray[2][0] +
		A.w * Other.pMatrixArray[3][0];

	Product.y = A.x * Other.pMatrixArray[0][1] +
		A.y * Other.pMatrixArray[1][1] +
		A.z * Other.pMatrixArray[2][1] +
		A.w * Other.pMatrixArray[3][1];

	Product.z = A.x * Other.pMatrixArray[0][2] +
		A.y * Other.pMatrixArray[1][2] +
		A.z * Other.pMatrixArray[2][2] +
		A.w * Other.pMatrixArray[3][2];

	Product.w = A.x * Other.pMatrixArray[0][3] +
		A.y * Other.pMatrixArray[1][3] +
		A.z * Other.pMatrixArray[2][3] +
		A.w * Other.pMatrixArray[3][3];
	return Product;
}


/*@note: Transforms Vertices by the projection and view matrix, the world matrix is not included as the identity matrix is used for the world
at the moment, not very ideal when multiple models ought to be drawn!
@param x: X origin
@param y: Y origin
@param nWidth: Viewport width
@param nHeight: Viewport height
@param View: View matrix
@param Projection: Projection matrix
@param Vertices: Vertices to be transformed*/
inline void To3D(int32_t x, int32_t y, int32_t nWidth, int32_t nHeight,
	MMatrix4x4& View, MMatrix4x4& Projection, MVector4* Vertices)
{
	MMatrix4x4 MVP = Projection * View;
	for (int i = 0; i < 3; i++)
		Vertices[i] = MVP * Vertices[i];
	for (int i = 0; i < 3; i++)
	{
		Vertices[i].x /= Vertices[i].w;
		Vertices[i].y /= Vertices[i].w;
		Vertices[i].z /= Vertices[i].w;
		Vertices[i].w /= Vertices[i].w;
	}
	for (int i = 0; i < 3; i++)
	{
		//@todo: fix this, passing same parameter is asking for trouble
		ViewportTransform(x, y, nWidth, nHeight, Vertices[i], Vertices[i]);
	}
}


namespace Math{

	/*@note: Performs a linear interpolation between A and B - specified by alpha*/
	template<typename T> T Math::Lerp(T& A, T& B, float Alpha)
	{
		return A * (1 - Alpha) + B * Alpha;
	}
	/*@note: Returns the absolute value of Type*/
	template<typename T> T Abs(T Type);
};