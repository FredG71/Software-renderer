/**
/*
*	Author: Frederic Garnier.
*	Language: C++
*	Created: 22/10/2014
*	Modified: 24/10/2014
*
*	
*
*	Notes:
*	Due to the way matrices and vectors were handled previously, a new matrix class and vector class
*	has been written, support for 2D vectors, utility functions for conversion between different coordinate systems
*	might be added as well. Matrices can now be inverted, transposed, cross product is implemented, functions to get the 
*	determinant has been implemented as well, some names have been changed for clarity. Multiplication between
*	a row or column vector is more explicit, operator overloading for multiplication has not and WILL not be implemented
*	to make things clearer, and avoid bugs. 
*	
*	Compiler: Microsoft Visual C++ Compiler Nov 2013 CTP
*
*	Platform tested: Windows 7 x64
*
*/

#include <math.h>
#include <cstdio>
#include <cfloat>
#include <smmintrin.h>
#include <intrin.h>
#include <cstdint>
#include <cassert>
#include <cstring>

#pragma warning( disable : 4201 )
typedef __m128 PackedVector;
typedef float Float4x4[4][4];
typedef float Float16[16];
typedef float FLOAT_PADDING;


#define SSE4_FLAGS 0x80000
#define ALIGN_16 __declspec( align( 16 ) )


#define LoadVector( X, Y, Z, W ) _mm_set_ps( W, Z, Y, X )
#define DotProduct( VectorA, VectorB, Mask ) _mm_dp_ps( VectorA, VectorB, Mask )
#define SetZero() _mm_setzero_ps()
#define SquareRoot( A ) _mm_sqrt_ps( A )
#define Divide( Dividend, Divisor ) _mm_div_ps(Dividend, Divisor )
#define Store( To, From ) _mm_store_ps( To, From );


/*Hints to be used for optimization*/
enum class EMatrix_Hints{
	// Matrix is a translation matrix
	ETranslation_Matrix,
	// Matrix is an orthogonal matrix
	EOrthogonal_Matrix,
	// Matrix is a scale matrix
	EScale_Matrix,
	// Unknown type - just get its inverse
	EUnknown_Matrix
};

enum class EMatrix_Translation_Properties{
	// Matrix transform(ed) a row vector
	ERow_Transform,
	// Matrix transform(ed) a column vector
	EColumn_Transform,
	// No translation in this matrix
	ENo_Translation
};

enum class EMatrix_Scale_Properties{
	// Scale transform is uniform
	EUniform_Scale,
	// Scale transform is non-uniorm
	ENonUniform_Scale,
	// No scale in this matrix
	ENo_Scale
};

namespace SSE{
	union UPackedMatrix{
		float fPacked[4][4];
		ALIGN_16 struct{
			PackedVector Row1;
			PackedVector Row2;
			PackedVector Row3;
			PackedVector Row4;
		};
	};
	union UPackedVector{
		ALIGN_16 struct{
			float fPacked[4];
		};
		PackedVector Vector;
	};
}

namespace MathUtility{
	/*@note: Flips the sign bit of fTarget - makes use of type punning,
	perhaps use a union or the reinterpret_cast operator to make things clear?
	@param fTarget: L-value reference to the float variable that'll have its sign flipped*/
	void Sign(float& fTarget)
	{
		(*(int32_t*)&fTarget) ^= 0x80000000;
	}
	bool IsEqual(float x, float y)
	{
		return abs(x - y) <= FLT_EPSILON * abs( x );
	}

	bool IsEqual(double x, double y)
	{
		return abs(x - y) <= DBL_EPSILON * abs( x );
	}
};


/*@note: Checks if SSE4 instructions are supported
@return: false if not supported, true otherwise*/
inline bool SSE4Support()
{
	int32_t pCPUInfo[4] { 0xFFFFFFFF };
	__cpuid(pCPUInfo, 0);
	if (pCPUInfo[0] < 1)
		return false;

	__cpuid(pCPUInfo, 1);

	return ((pCPUInfo[2] & SSE4_FLAGS) == SSE4_FLAGS);
}

class ALIGN_16 Vector4{
public:
	// Sets x, y, z, w to zero
	Vector4()
		: xmmVector(SetZero())
	{}

	// Sets x, y, z, w to the arguments passed w is 1 by default
	Vector4(float x, float y, float z, float w = 1)
		: xmmVector(LoadVector(x, y, z, w))
	{}

	// Normalizes this vector
	inline Vector4& Normalize()
	{
		// Get the square root of the dot product - 4 component
		PackedVector mDot = SquareRoot(DotProduct(xmmVector, xmmVector, 0xFF));
		// divide returned result
		xmmVector = Divide(xmmVector, mDot);
		return *this;
	}
	// returns the magnitude (length) of this vector
	inline float Length()
	{
		ALIGN_16 float fRetVal[4];
		// Gets the square root of the dot product
		PackedVector mDot = SquareRoot(DotProduct(xmmVector, xmmVector, 0xFF));
		// Store it and return it
		Store(fRetVal, mDot);
		return fRetVal[0];
	}
	union{
		ALIGN_16 struct{
			float x;
			float y;
			float z;
			float w;
		};
		PackedVector xmmVector;
	};
};

/* 3 Dimensional vector, fPadding is reserved*/
class ALIGN_16 Vector3{
public:
	Vector3()
		: xmmVector(SetZero())
	{}

	Vector3(float x, float y, float z)
		: xmmVector(LoadVector(x, y, z, 0))
	{}

	inline Vector3& Normalize()
	{
		PackedVector mDot = SquareRoot(DotProduct(xmmVector, xmmVector, 0x7F));
		xmmVector = Divide(xmmVector, mDot);
		return *this;
	}

	inline float Length()
	{
		ALIGN_16 float fRetVal[4];
		PackedVector mDot = SquareRoot(DotProduct(xmmVector, xmmVector, 0x7F));
		Store(fRetVal, mDot);
		return fRetVal[0];
	}
	union{
		ALIGN_16 struct{
			float x;
			float y;
			float z;
			FLOAT_PADDING __fPadding;
		};
		PackedVector xmmVector;
	};
};


class ALIGN_16 Matrix4x4{
public:
	explicit Matrix4x4(float fValue)
		: Row1(LoadVector(fValue, fValue, fValue, fValue)),
		Row2(LoadVector(fValue, fValue, fValue, fValue)),
		Row3(LoadVector(fValue, fValue, fValue, fValue)),
		Row4(LoadVector(fValue, fValue, fValue, fValue))
	{}
	Matrix4x4()
		: Row1(LoadVector(1, 0, 0, 0)),
		Row2(LoadVector(0, 1, 0, 0)),
		Row3(LoadVector(0, 0, 1, 0)),
		Row4(LoadVector(0, 0, 0, 1))
	{}
	Matrix4x4(float fRow1Column1, float fRow1Column2, float fRow1Column3, float fRow1Column4,
		float fRow2Column1, float fRow2Column2, float fRow2Column3, float fRow2Column4,
		float fRow3Column1, float fRow3Column2, float fRow3Column3, float fRow3Column4,
		float fRow4Column1, float fRow4Column2, float fRow4Column3, float fRow4Column4
		)
		: Row1(LoadVector(fRow1Column1, fRow1Column2, fRow1Column3, fRow1Column4)),
		Row2(LoadVector(fRow2Column1, fRow2Column2, fRow2Column3, fRow2Column4)),
		Row3(LoadVector(fRow3Column1, fRow3Column2, fRow3Column3, fRow3Column4)),
		Row4(LoadVector(fRow4Column1, fRow4Column2, fRow4Column3, fRow4Column4))
	{}

	inline Matrix4x4& operator *(float fScalar)
	{
		for (uint32_t nRow = 0; nRow <= 3; nRow++)
		{
			for (uint32_t nColumn = 0; nColumn <= 3; nColumn++)
			{
				MatrixArray[nRow][nColumn] *= fScalar;
			}
		}
		return *this;
	}
	/*@note: Gets a matrix component for the row nRow and the column nColumn, performs bounds checking
	@param nRow: Row of the matrix
	@param nColumn: Column of the matrix
	@return: Component in row nRow and column nColumn*/
	__forceinline float Get( uint32_t nRow, uint32_t nColumn )
	{
		assert( nRow <= 3 && nRow >= 0 && nColumn <= 3 && nColumn >= 0);
		return MatrixArray[nRow][nColumn];
	}

	/* @note: Computes the inverse of this matrix and gathers a few hints from the caller
		in cases where the inversion can be faster, for instance orthogonal matrices etc
	@param Hints: Hints to speed up the inversion
	@param TranslationProperties: What kind of translation is occuring, is it transforming
	a row or column vector?
	@param ScaleProperties: Is it a uniform or non-uniform scale?
	@return: Inversed matrix*/
	Matrix4x4& Inverse(
		EMatrix_Hints Hints, 
		EMatrix_Translation_Properties TranslationProperties = EMatrix_Translation_Properties::ERow_Transform,
		EMatrix_Scale_Properties ScaleProperties = EMatrix_Scale_Properties::EUniform_Scale
	)
	{
		switch (Hints)
		{
		case EMatrix_Hints::EOrthogonal_Matrix:
			return this->Transpose();
		case EMatrix_Hints::ETranslation_Matrix:
			switch (TranslationProperties)
			{
			case EMatrix_Translation_Properties::EColumn_Transform:
				MathUtility::Sign(MatrixArray[0][3]);
				MathUtility::Sign(MatrixArray[1][3]);
				MathUtility::Sign(MatrixArray[2][3]);
				break;
			case EMatrix_Translation_Properties::ERow_Transform:
				MathUtility::Sign(MatrixArray[3][1]);
				MathUtility::Sign(MatrixArray[3][2]);
				MathUtility::Sign(MatrixArray[3][3]);
				break;
			case EMatrix_Translation_Properties::ENo_Translation:
				break;
			}
			break;
		case EMatrix_Hints::EScale_Matrix:
			switch (ScaleProperties)
			{
				case EMatrix_Scale_Properties::EUniform_Scale:
				{
					register float fInverse =  1 / MatrixArray[0][0];
					MatrixArray[0][0] = fInverse;
					MatrixArray[1][1] = fInverse;
					MatrixArray[2][2] = fInverse;
					MatrixArray[3][3] = fInverse;
					return *this;
				}
				case EMatrix_Scale_Properties::ENonUniform_Scale:
				{
					register float fInverse = 1 / MatrixArray[0][0];
					MatrixArray[0][0] = fInverse;

					fInverse = 1 / MatrixArray[1][1];
					MatrixArray[1][1] = fInverse;

					fInverse = 1 / MatrixArray[2][2];
					MatrixArray[2][2] = fInverse;

					fInverse = 1 / MatrixArray[3][3];
					MatrixArray[3][3] = fInverse;
					return *this;
				}
			case EMatrix_Scale_Properties::ENo_Scale:
				break;
			}
			break;
		case EMatrix_Hints::EUnknown_Matrix:
				Float4x4 Result;
				Float4x4 Minors;
				float pDeterminant[4];
				Minors[0][0] = MatrixArray[2][2] * MatrixArray[3][3] - MatrixArray[2][3] * MatrixArray[3][2];
				Minors[0][1] = MatrixArray[1][2] * MatrixArray[3][3] - MatrixArray[1][3] * MatrixArray[3][2];
				Minors[0][2] = MatrixArray[1][2] * MatrixArray[2][3] - MatrixArray[1][3] * MatrixArray[2][2];

				Minors[1][0] = MatrixArray[2][2] * MatrixArray[3][3] - MatrixArray[2][3] * MatrixArray[3][2];
				Minors[1][1] = MatrixArray[0][2] * MatrixArray[3][3] - MatrixArray[0][3] * MatrixArray[3][2];
				Minors[1][2] = MatrixArray[0][2] * MatrixArray[2][3] - MatrixArray[0][3] * MatrixArray[2][2];

				Minors[2][0] = MatrixArray[1][2] * MatrixArray[3][3] - MatrixArray[1][3] * MatrixArray[3][2];
				Minors[2][1] = MatrixArray[0][2] * MatrixArray[3][3] - MatrixArray[0][3] * MatrixArray[3][2];
				Minors[2][2] = MatrixArray[0][2] * MatrixArray[1][3] - MatrixArray[0][3] * MatrixArray[1][2];

				Minors[3][0] = MatrixArray[1][2] * MatrixArray[2][3] - MatrixArray[1][3] * MatrixArray[2][2];
				Minors[3][1] = MatrixArray[0][2] * MatrixArray[2][3] - MatrixArray[0][3] * MatrixArray[2][2];
				Minors[3][2] = MatrixArray[0][2] * MatrixArray[1][3] - MatrixArray[0][3] * MatrixArray[1][2];

				pDeterminant[0] = MatrixArray[1][1] * Minors[0][0] - MatrixArray[2][1] * Minors[0][1] + MatrixArray[3][1] * Minors[0][2];
				pDeterminant[1] = MatrixArray[0][1] * Minors[1][0] - MatrixArray[2][1] * Minors[1][1] + MatrixArray[3][1] * Minors[1][2];
				pDeterminant[2] = MatrixArray[0][1] * Minors[2][0] - MatrixArray[1][1] * Minors[2][1] + MatrixArray[3][1] * Minors[2][2];
				pDeterminant[3] = MatrixArray[0][1] * Minors[3][0] - MatrixArray[1][1] * Minors[3][1] + MatrixArray[2][1] * Minors[3][2];

				float Determinant = MatrixArray[0][0] * pDeterminant[0] - MatrixArray[1][0] * pDeterminant[1] + MatrixArray[2][0] * pDeterminant[2] - MatrixArray[3][0] * pDeterminant[3];
				assert(!(MathUtility::IsEqual(Determinant, 0)));
				const float	RDet = 1.0f / Determinant;

				Result[0][0] = RDet * pDeterminant[0];
				Result[0][1] = -RDet * pDeterminant[1];
				Result[0][2] = RDet * pDeterminant[2];
				Result[0][3] = -RDet * pDeterminant[3];
				Result[1][0] = -RDet * (MatrixArray[1][0] * Minors[0][0] - MatrixArray[2][0] * Minors[0][1] + MatrixArray[3][0] * Minors[0][2]);
				Result[1][1] = RDet * (MatrixArray[0][0] * Minors[1][0] - MatrixArray[2][0] * Minors[1][1] + MatrixArray[3][0] * Minors[1][2]);
				Result[1][2] = -RDet * (MatrixArray[0][0] * Minors[2][0] - MatrixArray[1][0] * Minors[2][1] + MatrixArray[3][0] * Minors[2][2]);
				Result[1][3] = RDet * (MatrixArray[0][0] * Minors[3][0] - MatrixArray[1][0] * Minors[3][1] + MatrixArray[2][0] * Minors[3][2]);
				Result[2][0] = RDet * (
					MatrixArray[1][0] * (MatrixArray[2][1] * MatrixArray[3][3] - MatrixArray[2][3] * MatrixArray[3][1]) -
					MatrixArray[2][0] * (MatrixArray[1][1] * MatrixArray[3][3] - MatrixArray[1][3] * MatrixArray[3][1]) +
					MatrixArray[3][0] * (MatrixArray[1][1] * MatrixArray[2][3] - MatrixArray[1][3] * MatrixArray[2][1])
					);
				Result[2][1] = -RDet * (
					MatrixArray[0][0] * (MatrixArray[2][1] * MatrixArray[3][3] - MatrixArray[2][3] * MatrixArray[3][1]) -
					MatrixArray[2][0] * (MatrixArray[0][1] * MatrixArray[3][3] - MatrixArray[0][3] * MatrixArray[3][1]) +
					MatrixArray[3][0] * (MatrixArray[0][1] * MatrixArray[2][3] - MatrixArray[0][3] * MatrixArray[2][1])
					);
				Result[2][2] = RDet * (
					MatrixArray[0][0] * (MatrixArray[1][1] * MatrixArray[3][3] - MatrixArray[1][3] * MatrixArray[3][1]) -
					MatrixArray[1][0] * (MatrixArray[0][1] * MatrixArray[3][3] - MatrixArray[0][3] * MatrixArray[3][1]) +
					MatrixArray[3][0] * (MatrixArray[0][1] * MatrixArray[1][3] - MatrixArray[0][3] * MatrixArray[1][1])
					);
				Result[2][3] = -RDet * (
					MatrixArray[0][0] * (MatrixArray[1][1] * MatrixArray[2][3] - MatrixArray[1][3] * MatrixArray[2][1]) -
					MatrixArray[1][0] * (MatrixArray[0][1] * MatrixArray[2][3] - MatrixArray[0][3] * MatrixArray[2][1]) +
					MatrixArray[2][0] * (MatrixArray[0][1] * MatrixArray[1][3] - MatrixArray[0][3] * MatrixArray[1][1])
					);
				Result[3][0] = -RDet * (
					MatrixArray[1][0] * (MatrixArray[2][1] * MatrixArray[3][2] - MatrixArray[2][2] * MatrixArray[3][1]) -
					MatrixArray[2][0] * (MatrixArray[1][1] * MatrixArray[3][2] - MatrixArray[1][2] * MatrixArray[3][1]) +
					MatrixArray[3][0] * (MatrixArray[1][1] * MatrixArray[2][2] - MatrixArray[1][2] * MatrixArray[2][1])
					);
				Result[3][1] = RDet * (
					MatrixArray[0][0] * (MatrixArray[2][1] * MatrixArray[3][2] - MatrixArray[2][2] * MatrixArray[3][1]) -
					MatrixArray[2][0] * (MatrixArray[0][1] * MatrixArray[3][2] - MatrixArray[0][2] * MatrixArray[3][1]) +
					MatrixArray[3][0] * (MatrixArray[0][1] * MatrixArray[2][2] - MatrixArray[0][2] * MatrixArray[2][1])
					);
				Result[3][2] = -RDet * (
					MatrixArray[0][0] * (MatrixArray[1][1] * MatrixArray[3][2] - MatrixArray[1][2] * MatrixArray[3][1]) -
					MatrixArray[1][0] * (MatrixArray[0][1] * MatrixArray[3][2] - MatrixArray[0][2] * MatrixArray[3][1]) +
					MatrixArray[3][0] * (MatrixArray[0][1] * MatrixArray[1][2] - MatrixArray[0][2] * MatrixArray[1][1])
					);
				Result[3][3] = RDet * (
					MatrixArray[0][0] * (MatrixArray[1][1] * MatrixArray[2][2] - MatrixArray[1][2] * MatrixArray[2][1]) -
					MatrixArray[1][0] * (MatrixArray[0][1] * MatrixArray[2][2] - MatrixArray[0][2] * MatrixArray[2][1]) +
					MatrixArray[2][0] * (MatrixArray[0][1] * MatrixArray[1][2] - MatrixArray[0][2] * MatrixArray[1][1])
					);

				memcpy(MatrixArray, &Result, sizeof(Float4x4));
			break;
		}
		return *this;
	}
	/*@note: Returns the determiment of this matrix
	@param: Matrix determinent*/
	inline float GetDeterminant()
	{
		return	MatrixArray[0][0] * (
				MatrixArray[1][1] * (MatrixArray[2][2] * MatrixArray[3][3] - MatrixArray[2][3] * MatrixArray[3][2]) -
				MatrixArray[2][1] * (MatrixArray[1][2] * MatrixArray[3][3] - MatrixArray[1][3] * MatrixArray[3][2]) +
				MatrixArray[3][1] * (MatrixArray[1][2] * MatrixArray[2][3] - MatrixArray[1][3] * MatrixArray[2][2]))

				-

				MatrixArray[1][0] * (
				MatrixArray[0][1] * (MatrixArray[2][2] * MatrixArray[3][3] - MatrixArray[2][3] * MatrixArray[3][2]) -
				MatrixArray[2][1] * (MatrixArray[0][2] * MatrixArray[3][3] - MatrixArray[0][3] * MatrixArray[3][2]) +
				MatrixArray[3][1] * (MatrixArray[0][2] * MatrixArray[2][3] - MatrixArray[0][3] * MatrixArray[2][2])) 
				
				+

				MatrixArray[2][0] * (
				MatrixArray[0][1] * (MatrixArray[1][2] * MatrixArray[3][3] - MatrixArray[1][3] * MatrixArray[3][2]) -
				MatrixArray[1][1] * (MatrixArray[0][2] * MatrixArray[3][3] - MatrixArray[0][3] * MatrixArray[3][2]) +
				MatrixArray[3][1] * (MatrixArray[0][2] * MatrixArray[1][3] - MatrixArray[0][3] * MatrixArray[1][2])) 
				
				-

				MatrixArray[3][0] * (
				MatrixArray[0][1] * (MatrixArray[1][2] * MatrixArray[2][3] - MatrixArray[1][3] * MatrixArray[2][2]) -
				MatrixArray[1][1] * (MatrixArray[0][2] * MatrixArray[2][3] - MatrixArray[0][3] * MatrixArray[2][2]) +
				MatrixArray[2][1] * (MatrixArray[0][2] * MatrixArray[1][3] - MatrixArray[0][3] * MatrixArray[1][2]));
	}
	/*@note: Transposes this matrix, note that (this) is modified
	@return: Transposed matrix*/
	inline Matrix4x4& Transpose()
	{
		// Load columns into PackedVectors
		PackedVector Column1 = LoadVector(Row1X, Row2X, Row3X, Row4X);
		PackedVector Column2 = LoadVector(Row1Y, Row2Y, Row3Y, Row4Y);
		PackedVector Column3 = LoadVector(Row1Z, Row2Z, Row3Z, Row4Z);
		PackedVector Column4 = LoadVector(Row1W, Row2W, Row3W, Row4W);

		// Switch rows with columns
		Row1 = Column1;
		Row2 = Column2;
		Row3 = Column3;
		Row4 = Column4;

		return *this;
	}
	/*@note: Transforms a column vector by this matrix. This is shorter due to the way
	vectors are stored.
	@param ColumnVector: R-value reference to the column vector that is to be transformed by this matrix
	@return: Transformed vector*/
	Vector4 TransformColumnVector(Vector4&& ColumnVector)
	{
		Matrix4x4 Temp = {};
		Vector4 RetVal = {};

		// Get the dot product of the rows of Temp (this) with the columns of ColumnVector
		Temp.Row1 = DotProduct(Row1, ColumnVector.xmmVector, 0xFF);
		Temp.Row2 = DotProduct(Row2, ColumnVector.xmmVector, 0xFF);
		Temp.Row3 = DotProduct(Row3, ColumnVector.xmmVector, 0xFF);
		Temp.Row4 = DotProduct(Row4, ColumnVector.xmmVector, 0xFF);

		// Get a packed union
		SSE::UPackedVector PackedUnion = {};
		// And cache the results to a PackedVector
		PackedUnion.Vector = Temp.Row1;
		RetVal.x = PackedUnion.fPacked[0];

		PackedUnion.Vector = Temp.Row2;
		RetVal.y = PackedUnion.fPacked[1];

		PackedUnion.Vector = Temp.Row3;
		RetVal.z = PackedUnion.fPacked[2];

		PackedUnion.Vector = Temp.Row4;
		RetVal.w = PackedUnion.fPacked[3];

		return RetVal;
	}
	/*@note: Transforms a column vector by this matrix. This is shorter due to the way
	vectors are stored.
	@param ColumnVector: L-value reference to the column vector that is to be transformed by this matrix
	@return: Transformed vector*/
	Vector4 TransformColumnVector(Vector4& ColumnVector)
	{
		Matrix4x4 Temp = {};
		Vector4 RetVal = {};

		Temp.Row1 = DotProduct(Row1, ColumnVector.xmmVector, 0xFF);
		Temp.Row2 = DotProduct(Row2, ColumnVector.xmmVector, 0xFF);
		Temp.Row3 = DotProduct(Row3, ColumnVector.xmmVector, 0xFF);
		Temp.Row4 = DotProduct(Row4, ColumnVector.xmmVector, 0xFF);

		SSE::UPackedVector PackedUnion = {};
		PackedUnion.Vector = Temp.Row1;
		RetVal.x = PackedUnion.fPacked[0];

		PackedUnion.Vector = Temp.Row2;
		RetVal.y = PackedUnion.fPacked[1];

		PackedUnion.Vector = Temp.Row3;
		RetVal.z = PackedUnion.fPacked[2];

		PackedUnion.Vector = Temp.Row4;
		RetVal.w = PackedUnion.fPacked[3];

		return RetVal;
	}
	/*@note: Transforms a row vector by this matrix
	@param Row vector: R-value reference to the row vector that is to be transformed by this matrix
	@return: Transformed vector*/
	Vector4 TransformRowVector(Vector4&& RowVector)
	{

		Vector4 RetVal = {};
		// Get the columns of this matrix and store them in a PackedVector
		PackedVector Column1 = LoadVector(Row1X, Row2X, Row3X, Row4X);
		PackedVector Column2 = LoadVector(Row1Y, Row2Y, Row3Y, Row4Y);
		PackedVector Column3 = LoadVector(Row1Z, Row2Z, Row3Z, Row4Z);
		PackedVector Column4 = LoadVector(Row1W, Row2W, Row3W, Row4W);


		// Get the dot product of the rows of RowVector and the columns of this matrix
		PackedVector Row1Column1 = DotProduct(RowVector.xmmVector, Column1, 0xFF);
		PackedVector Row1Column2 = DotProduct(RowVector.xmmVector, Column2, 0xFF);
		PackedVector Row1Column3 = DotProduct(RowVector.xmmVector, Column3, 0xFF);
		PackedVector Row1Column4 = DotProduct(RowVector.xmmVector, Column4, 0xFF);

		// Get a packed union
		SSE::UPackedVector PackedUnion = {};

		// And store the results accordingly
		PackedUnion.Vector = Row1Column1;
		RetVal.x = PackedUnion.fPacked[0];

		PackedUnion.Vector = Row1Column2;
		RetVal.y = PackedUnion.fPacked[1];

		PackedUnion.Vector = Row1Column3;
		RetVal.z = PackedUnion.fPacked[2];

		PackedUnion.Vector = Row1Column4;
		RetVal.w = PackedUnion.fPacked[3];

		return RetVal;
	}
	/*@note: Transforms a row vector by this matrix
	@param Row vector: L-value reference to the row vector that is to be transformed by this matrix
	@return: Transformed vector*/
	Vector4 TransformRowVector(Vector4& RowVector)
	{

		Vector4 RetVal = {};
		PackedVector Column1 = LoadVector(Row1X, Row2X, Row3X, Row4X);
		PackedVector Column2 = LoadVector(Row1Y, Row2Y, Row3Y, Row4Y);
		PackedVector Column3 = LoadVector(Row1Z, Row2Z, Row3Z, Row4Z);
		PackedVector Column4 = LoadVector(Row1W, Row2W, Row3W, Row4W);

		PackedVector Row1Column1 = DotProduct(RowVector.xmmVector, Column1, 0xFF);
		PackedVector Row1Column2 = DotProduct(RowVector.xmmVector, Column2, 0xFF);
		PackedVector Row1Column3 = DotProduct(RowVector.xmmVector, Column3, 0xFF);
		PackedVector Row1Column4 = DotProduct(RowVector.xmmVector, Column4, 0xFF);

		SSE::UPackedVector PackedUnion = {};

		PackedUnion.Vector = Row1Column1;
		RetVal.x = PackedUnion.fPacked[0];

		PackedUnion.Vector = Row1Column2;
		RetVal.y = PackedUnion.fPacked[1];

		PackedUnion.Vector = Row1Column3;
		RetVal.z = PackedUnion.fPacked[2];

		PackedUnion.Vector = Row1Column4;
		RetVal.w = PackedUnion.fPacked[3];

		return RetVal;
	}
	union{
		ALIGN_16 struct{
			float Row1X;
			float Row1Y;
			float Row1Z;
			float Row1W;

			float Row2X;
			float Row2Y;
			float Row2Z;
			float Row2W;

			float Row3X;
			float Row3Y;
			float Row3Z;
			float Row3W;

			float Row4X;
			float Row4Y;
			float Row4Z;
			float Row4W;
		};
		Float4x4 MatrixArray;
		ALIGN_16 struct{
			PackedVector Row1;
			PackedVector Row2;
			PackedVector Row3;
			PackedVector Row4;
		};
	};
};

namespace Math{

	/*@note: Returns the length (magnitude) of A
	@param A: R-value reference to a Vector4
	@return: The magnitude of A*/
	inline float Length(Vector4&& A)
	{
		ALIGN_16 float fRetVal[4];
		// Get the magnitude
		PackedVector mDot = SquareRoot(DotProduct(A.xmmVector, A.xmmVector, 0xFF));
		// Store it
		Store(fRetVal, mDot);
		// And return it
		return fRetVal[0];
	}
	/*@note: Returns the length (magnitude) of A
	@param A: L-value reference to a Vector4
	@return: The magnitude of A*/
	inline float Length(Vector4& A)
	{
		ALIGN_16 float fRetVal[4];
		PackedVector mDot = SquareRoot(DotProduct(A.xmmVector, A.xmmVector, 0xFF));
		Store(fRetVal, mDot);
		return fRetVal[0];
	}

	/*@note: Returns the length (magnitude) of A
	@param A: R-value reference to a Vector3
	@return: The magnitude of A*/
	inline float Length(Vector3&& A)
	{
		ALIGN_16 float fRetVal[4];
		// Get the magnitude - ensure no division by zero occurs
		PackedVector mDot = SquareRoot(DotProduct(A.xmmVector, A.xmmVector, 0x7F));
		// Store result
		Store(fRetVal, mDot);
		// And return it
		return fRetVal[0];
	}

	/*@note: Returns the length (magnitude) of A
	@param A: L-value reference to a Vector3
	@return: The magnitude of A*/
	inline float Length(Vector3& A)
	{
		ALIGN_16 float fRetVal[4];
		PackedVector mDot = SquareRoot(DotProduct(A.xmmVector, A.xmmVector, 0x7F));
		Store(fRetVal, mDot);
		return fRetVal[0];
	}

	/*note: Returns the dot product of A and B
	@param A: R-value reference to a Vector4
	@param B: R-Value reference to a Vector4
	@return: The resulted dot product*/
	inline float Dot(Vector4&& A, Vector4&& B)
	{
		ALIGN_16 float fResult[4];
		PackedVector RetVal = DotProduct(A.xmmVector, B.xmmVector, 0xFF); // 1111 1111 - sum 4 components
		Store(fResult, RetVal);
		return fResult[0];
	}
	/*note: Returns the dot product of A and B
	@param A: L-value reference to a Vector4
	@param B: R-Value reference to a Vector4
	@return: The resulted dot product*/
	inline float Dot(Vector4& A, Vector4&& B)
	{
		ALIGN_16 float fResult[4];
		PackedVector RetVal = DotProduct(A.xmmVector, B.xmmVector, 0xFF); // 1111 1111 - sum 4 components
		Store(fResult, RetVal);
		return fResult[0];
	}
	/*note: Returns the dot product of A and B
	@param A: R-value reference to a Vector4
	@param B: L-Value reference to a Vector4
	@return: The resulted dot product*/
	inline float Dot(Vector4&& A, Vector4& B)
	{
		ALIGN_16 float fResult[4];
		PackedVector RetVal = DotProduct(A.xmmVector, B.xmmVector, 0xFF); // 1111 1111 - sum 4 components
		Store(fResult, RetVal);
		return fResult[0];
	}
	/*note: Returns the dot product of A and B
	@param A: L-value reference to a Vector4
	@param B: L-Value reference to a Vector4
	@return: The resulted dot product*/
	inline float Dot(Vector4& A, Vector4& B)
	{
		ALIGN_16 float fResult[4];
		PackedVector RetVal = DotProduct(A.xmmVector, B.xmmVector, 0xFF); // 1111 1111 - sum 4 components
		Store(fResult, RetVal);
		return fResult[0];
	}

	/*note: Returns the dot product of A and B
	@param A: R-value reference to a Vector3
	@param B: R-Value reference to a Vector3
	@return: The resulted dot product*/
	inline float Dot(Vector3&& A, Vector3&& B)
	{
		ALIGN_16 float fResult[4];
		PackedVector RetVal = DotProduct(A.xmmVector, B.xmmVector, 0x71); // 0111 0001 - only sum 3 components
		Store(fResult, RetVal);
		return fResult[0];
	}
	/*note: Returns the dot product of A and B
	@param A: L-value reference to a Vector3
	@param B: R-Value reference to a Vector3
	@return: The resulted dot product*/
	inline float Dot(Vector3& A, Vector3&& B)
	{
		ALIGN_16 float fResult[4];
		PackedVector RetVal = DotProduct(A.xmmVector, B.xmmVector, 0x71); // 0111 0001 - only sum 3 components
		Store(fResult, RetVal);
		return fResult[0];
	}
	/*note: Returns the dot product of A and B
	@param A: R-value reference to a Vector3
	@param B: L-Value reference to a Vector3
	@return: The resulted dot product*/
	inline float Dot(Vector3&& A, Vector3& B)
	{
		ALIGN_16 float fResult[4];
		PackedVector RetVal = DotProduct(A.xmmVector, B.xmmVector, 0x71); // 0111 0001 - only sum 3 components
		Store(fResult, RetVal);
		return fResult[0];
	}
	/*note: Returns the dot product of A and B
	@param A: L-value reference to a Vector3
	@param B: L-Value reference to a Vector3
	@return: The resulted dot product*/
	inline float Dot(Vector3& A, Vector3& B)
	{
		ALIGN_16 float fResult[4];
		PackedVector RetVal = DotProduct(A.xmmVector, B.xmmVector, 0x71); // 0111 0001 - only sum 3 components
		Store(fResult, RetVal);
		return fResult[0];
	}

	/*@note: Returns the cross product of A and B - order A Cross B
	x86 assembly was used to make things short
	@param A: L-value reference - Left operand of the cross product
	@param B: L-value reference - Right operand of the cross product
	@return: Vector orthogonal to A and B*/
	inline Vector3 Cross(Vector3& A, Vector3& B)
	{
		Vector3 vRetVal = Vector3();
		__asm{
			mov eax, A
				mov ecx, B

				movaps xmm0, [eax]
				movaps xmm1, [ecx]

				movaps xmm2, xmm0
				movaps xmm3, xmm1; keep copies of A and B

				shufps xmm0, xmm0, 0xD8; 1101 1000b
				shufps xmm1, xmm1, 0xE1; 1110 0001b
				mulps xmm0, xmm1

				shufps xmm2, xmm2, 0xE1; 1110 0001b
				shufps xmm3, xmm3, 0xD8; 1101 1000b
				mulps xmm2, xmm3

				subps xmm0, xmm2
				shufps xmm0, xmm0, 0xC6; 1100 0110b (-not part of cross product switches x and z)
				movaps[vRetVal], xmm0
		}
		return vRetVal;
	}
	/*@note: Returns the cross product of A and B - order A Cross B
	x86 assembly was used to make things short
	@param A: R-value reference - Left operand of the cross product
	@param B: R-value reference - Right operand of the cross product
	@return: Vector orthogonal to A and B*/
	inline Vector3 Cross(Vector3&& A, Vector3&& B)
	{
		Vector3 vRetVal = Vector3();
		__asm{
			mov eax, A
				mov ecx, B

				movaps xmm0, [eax]
				movaps xmm1, [ecx]

				movaps xmm2, xmm0
				movaps xmm3, xmm1; keep copies of A and B

				shufps xmm0, xmm0, 0xD8; 1101 1000b
				shufps xmm1, xmm1, 0xE1; 1110 0001b
				mulps xmm0, xmm1

				shufps xmm2, xmm2, 0xE1; 1110 0001b
				shufps xmm3, xmm3, 0xD8; 1101 1000b
				mulps xmm2, xmm3

				subps xmm0, xmm2
				shufps xmm0, xmm0, 0xC6; 1100 0110b (-not part of cross product switches x and z)
				movaps[vRetVal], xmm0
		}
		return vRetVal;
	}
	/*@note: Returns the cross product of A and B - order A Cross B
	x86 assembly was used to make things short
	@param A: L-value reference - Left operand of the cross product
	@param B: R-value reference - Right operand of the cross product
	@return: Vector orthogonal to A and B*/
	inline Vector3 Cross(Vector3& A, Vector3&& B)
	{
		Vector3 vRetVal = Vector3();
		__asm{
			mov eax, A
				mov ecx, B

				movaps xmm0, [eax]
				movaps xmm1, [ecx]

				movaps xmm2, xmm0
				movaps xmm3, xmm1; keep copies of A and B

				shufps xmm0, xmm0, 0xD8; 1101 1000b
				shufps xmm1, xmm1, 0xE1; 1110 0001b
				mulps xmm0, xmm1

				shufps xmm2, xmm2, 0xE1; 1110 0001b
				shufps xmm3, xmm3, 0xD8; 1101 1000b
				mulps xmm2, xmm3

				subps xmm0, xmm2
				shufps xmm0, xmm0, 0xC6; 1100 0110b (-not part of cross product switches x and z)
				movaps[vRetVal], xmm0
		}
		return vRetVal;
	}
	/*@note: Returns the cross product of A and B - order A Cross B
	x86 assembly was used to make things short
	@param A: R-value reference - Left operand of the cross product
	@param B: L-value reference - Right operand of the cross product
	@return: Vector orthogonal to A and B*/
	inline Vector3 Cross(Vector3&& A, Vector3& B)
	{
		Vector3 vRetVal = Vector3();
		__asm{
			mov eax, A
				mov ecx, B

				movaps xmm0, [eax]
				movaps xmm1, [ecx]

				movaps xmm2, xmm0
				movaps xmm3, xmm1; keep copies of A and B

				shufps xmm0, xmm0, 0xD8; 1101 1000b
				shufps xmm1, xmm1, 0xE1; 1110 0001b
				mulps xmm0, xmm1

				shufps xmm2, xmm2, 0xE1; 1110 0001b
				shufps xmm3, xmm3, 0xD8; 1101 1000b
				mulps xmm2, xmm3

				subps xmm0, xmm2
				shufps xmm0, xmm0, 0xC6; 1100 0110b (-not part of cross product switches x and z)
				movaps[vRetVal], xmm0
		}
		return vRetVal;
	}


	/*@note: Multiplies Matrix A by Matrix B
	@param A: L-value reference to a 4 by 4 matrix
	@param B: L-value reference to a 4 by 4 matrix
	@return: Resulting product*/
	Matrix4x4 MatrixProduct(Matrix4x4& A, Matrix4x4& B)
	{
		Matrix4x4 Result = {};

		// Get the columns
		PackedVector Column1 = LoadVector(B.Row1X, B.Row2X, B.Row3X, B.Row4X);
		PackedVector Column2 = LoadVector(B.Row1Y, B.Row2Y, B.Row3Y, B.Row4Y);
		PackedVector Column3 = LoadVector(B.Row1Z, B.Row2Z, B.Row3Z, B.Row4Z);
		PackedVector Column4 = LoadVector(B.Row1W, B.Row2W, B.Row3W, B.Row4W);

		// And store the dot product of the rows of A with
		// the columns of B to PackedVectors

		/*First row*/
		PackedVector Row1Column1 = DotProduct(A.Row1, Column1, 0xFF);
		PackedVector Row1Column2 = DotProduct(A.Row1, Column2, 0xFF);
		PackedVector Row1Column3 = DotProduct(A.Row1, Column3, 0xFF);
		PackedVector Row1Column4 = DotProduct(A.Row1, Column4, 0xFF);

		/*Second row*/
		PackedVector Row2Column1 = DotProduct(A.Row2, Column1, 0xFF);
		PackedVector Row2Column2 = DotProduct(A.Row2, Column2, 0xFF);
		PackedVector Row2Column3 = DotProduct(A.Row2, Column3, 0xFF);
		PackedVector Row2Column4 = DotProduct(A.Row2, Column4, 0xFF);

		/*Third row*/
		PackedVector Row3Column1 = DotProduct(A.Row3, Column1, 0xFF);
		PackedVector Row3Column2 = DotProduct(A.Row3, Column2, 0xFF);
		PackedVector Row3Column3 = DotProduct(A.Row3, Column3, 0xFF);
		PackedVector Row3Column4 = DotProduct(A.Row3, Column4, 0xFF);

		/*Fourth row*/
		PackedVector Row4Column1 = DotProduct(A.Row4, Column1, 0xFF);
		PackedVector Row4Column2 = DotProduct(A.Row4, Column2, 0xFF);
		PackedVector Row4Column3 = DotProduct(A.Row4, Column3, 0xFF);
		PackedVector Row4Column4 = DotProduct(A.Row4, Column4, 0xFF);

		// Get a packed matrix - union
		SSE::UPackedMatrix PackedResult = {};

		// And store them accordingly - indexing isn't needed but is used to make things clear
		// First row
		PackedResult.Row1 = Row1Column1;
		Result.Row1X = PackedResult.fPacked[0][0];
		PackedResult.Row1 = Row1Column2;
		Result.Row1Y = PackedResult.fPacked[0][1];
		PackedResult.Row1 = Row1Column3;
		Result.Row1Z = PackedResult.fPacked[0][2];
		PackedResult.Row1 = Row1Column4;
		Result.Row1W = PackedResult.fPacked[0][3];

		// Second row
		PackedResult.Row2 = Row2Column1;
		Result.Row2X = PackedResult.fPacked[1][0];
		PackedResult.Row2 = Row2Column2;
		Result.Row2Y = PackedResult.fPacked[1][1];
		PackedResult.Row2 = Row2Column3;
		Result.Row2Z = PackedResult.fPacked[1][2];
		PackedResult.Row2 = Row2Column4;
		Result.Row2W = PackedResult.fPacked[1][3];

		// Third row
		PackedResult.Row3 = Row3Column1;
		Result.Row3X = PackedResult.fPacked[2][0];
		PackedResult.Row3 = Row3Column2;
		Result.Row3Y = PackedResult.fPacked[2][1];
		PackedResult.Row3 = Row3Column3;
		Result.Row3Z = PackedResult.fPacked[2][2];
		PackedResult.Row3 = Row3Column4;
		Result.Row3W = PackedResult.fPacked[2][3];

		PackedResult.Row4 = Row4Column1;
		Result.Row4X = PackedResult.fPacked[3][0];
		PackedResult.Row4 = Row4Column2;
		Result.Row4Y = PackedResult.fPacked[3][1];
		PackedResult.Row4 = Row4Column3;
		Result.Row4Z = PackedResult.fPacked[3][2];
		PackedResult.Row4 = Row4Column4;
		Result.Row4W = PackedResult.fPacked[3][3];

		return Result;
	}

	/*@note: Multiplies Matrix A by Matrix B
	@param A: L-value reference to a 4 by 4 matrix
	@param B: L-value reference to a 4 by 4 matrix
	@return: Resulting product*/
	Matrix4x4 MatrixProduct(Matrix4x4&& A, Matrix4x4& B)
	{
		Matrix4x4 Result = {};

		// Get the columns
		PackedVector Column1 = LoadVector(B.Row1X, B.Row2X, B.Row3X, B.Row4X);
		PackedVector Column2 = LoadVector(B.Row1Y, B.Row2Y, B.Row3Y, B.Row4Y);
		PackedVector Column3 = LoadVector(B.Row1Z, B.Row2Z, B.Row3Z, B.Row4Z);
		PackedVector Column4 = LoadVector(B.Row1W, B.Row2W, B.Row3W, B.Row4W);

		// And store the dot product of the rows of A with
		// the columns of B to PackedVectors

		/*First row*/
		PackedVector Row1Column1 = DotProduct(A.Row1, Column1, 0xFF);
		PackedVector Row1Column2 = DotProduct(A.Row1, Column2, 0xFF);
		PackedVector Row1Column3 = DotProduct(A.Row1, Column3, 0xFF);
		PackedVector Row1Column4 = DotProduct(A.Row1, Column4, 0xFF);

		/*Second row*/
		PackedVector Row2Column1 = DotProduct(A.Row2, Column1, 0xFF);
		PackedVector Row2Column2 = DotProduct(A.Row2, Column2, 0xFF);
		PackedVector Row2Column3 = DotProduct(A.Row2, Column3, 0xFF);
		PackedVector Row2Column4 = DotProduct(A.Row2, Column4, 0xFF);

		/*Third row*/
		PackedVector Row3Column1 = DotProduct(A.Row3, Column1, 0xFF);
		PackedVector Row3Column2 = DotProduct(A.Row3, Column2, 0xFF);
		PackedVector Row3Column3 = DotProduct(A.Row3, Column3, 0xFF);
		PackedVector Row3Column4 = DotProduct(A.Row3, Column4, 0xFF);

		/*Fourth row*/
		PackedVector Row4Column1 = DotProduct(A.Row4, Column1, 0xFF);
		PackedVector Row4Column2 = DotProduct(A.Row4, Column2, 0xFF);
		PackedVector Row4Column3 = DotProduct(A.Row4, Column3, 0xFF);
		PackedVector Row4Column4 = DotProduct(A.Row4, Column4, 0xFF);

		// Get a packed matrix - union
		SSE::UPackedMatrix PackedResult = {};

		// And store them accordingly - indexing isn't needed but is used to make things clear
		// First row
		PackedResult.Row1 = Row1Column1;
		Result.Row1X = PackedResult.fPacked[0][0];
		PackedResult.Row1 = Row1Column2;
		Result.Row1Y = PackedResult.fPacked[0][1];
		PackedResult.Row1 = Row1Column3;
		Result.Row1Z = PackedResult.fPacked[0][2];
		PackedResult.Row1 = Row1Column4;
		Result.Row1W = PackedResult.fPacked[0][3];

		// Second row
		PackedResult.Row2 = Row2Column1;
		Result.Row2X = PackedResult.fPacked[1][0];
		PackedResult.Row2 = Row2Column2;
		Result.Row2Y = PackedResult.fPacked[1][1];
		PackedResult.Row2 = Row2Column3;
		Result.Row2Z = PackedResult.fPacked[1][2];
		PackedResult.Row2 = Row2Column4;
		Result.Row2W = PackedResult.fPacked[1][3];

		// Third row
		PackedResult.Row3 = Row3Column1;
		Result.Row3X = PackedResult.fPacked[2][0];
		PackedResult.Row3 = Row3Column2;
		Result.Row3Y = PackedResult.fPacked[2][1];
		PackedResult.Row3 = Row3Column3;
		Result.Row3Z = PackedResult.fPacked[2][2];
		PackedResult.Row3 = Row3Column4;
		Result.Row3W = PackedResult.fPacked[2][3];

		PackedResult.Row4 = Row4Column1;
		Result.Row4X = PackedResult.fPacked[3][0];
		PackedResult.Row4 = Row4Column2;
		Result.Row4Y = PackedResult.fPacked[3][1];
		PackedResult.Row4 = Row4Column3;
		Result.Row4Z = PackedResult.fPacked[3][2];
		PackedResult.Row4 = Row4Column4;
		Result.Row4W = PackedResult.fPacked[3][3];

		return Result;
	}

	/*@note: Multiplies Matrix A by Matrix B
	@param A: L-value reference to a 4 by 4 matrix
	@param B: R-value reference to a 4 by 4 matrix
	@return: Resulting product*/
	Matrix4x4 MatrixProduct(Matrix4x4& A, Matrix4x4&& B)
	{
		Matrix4x4 Result = {};

		// Get the columns
		PackedVector Column1 = LoadVector(B.Row1X, B.Row2X, B.Row3X, B.Row4X);
		PackedVector Column2 = LoadVector(B.Row1Y, B.Row2Y, B.Row3Y, B.Row4Y);
		PackedVector Column3 = LoadVector(B.Row1Z, B.Row2Z, B.Row3Z, B.Row4Z);
		PackedVector Column4 = LoadVector(B.Row1W, B.Row2W, B.Row3W, B.Row4W);

		// And store the dot product of the rows of A with
		// the columns of B to PackedVectors

		/*First row*/
		PackedVector Row1Column1 = DotProduct(A.Row1, Column1, 0xFF);
		PackedVector Row1Column2 = DotProduct(A.Row1, Column2, 0xFF);
		PackedVector Row1Column3 = DotProduct(A.Row1, Column3, 0xFF);
		PackedVector Row1Column4 = DotProduct(A.Row1, Column4, 0xFF);

		/*Second row*/
		PackedVector Row2Column1 = DotProduct(A.Row2, Column1, 0xFF);
		PackedVector Row2Column2 = DotProduct(A.Row2, Column2, 0xFF);
		PackedVector Row2Column3 = DotProduct(A.Row2, Column3, 0xFF);
		PackedVector Row2Column4 = DotProduct(A.Row2, Column4, 0xFF);

		/*Third row*/
		PackedVector Row3Column1 = DotProduct(A.Row3, Column1, 0xFF);
		PackedVector Row3Column2 = DotProduct(A.Row3, Column2, 0xFF);
		PackedVector Row3Column3 = DotProduct(A.Row3, Column3, 0xFF);
		PackedVector Row3Column4 = DotProduct(A.Row3, Column4, 0xFF);

		/*Fourth row*/
		PackedVector Row4Column1 = DotProduct(A.Row4, Column1, 0xFF);
		PackedVector Row4Column2 = DotProduct(A.Row4, Column2, 0xFF);
		PackedVector Row4Column3 = DotProduct(A.Row4, Column3, 0xFF);
		PackedVector Row4Column4 = DotProduct(A.Row4, Column4, 0xFF);

		// Get a packed matrix - union
		SSE::UPackedMatrix PackedResult = {};

		// And store them accordingly - indexing isn't needed but is used to make things clear
		// First row
		PackedResult.Row1 = Row1Column1;
		Result.Row1X = PackedResult.fPacked[0][0];
		PackedResult.Row1 = Row1Column2;
		Result.Row1Y = PackedResult.fPacked[0][1];
		PackedResult.Row1 = Row1Column3;
		Result.Row1Z = PackedResult.fPacked[0][2];
		PackedResult.Row1 = Row1Column4;
		Result.Row1W = PackedResult.fPacked[0][3];

		// Second row
		PackedResult.Row2 = Row2Column1;
		Result.Row2X = PackedResult.fPacked[1][0];
		PackedResult.Row2 = Row2Column2;
		Result.Row2Y = PackedResult.fPacked[1][1];
		PackedResult.Row2 = Row2Column3;
		Result.Row2Z = PackedResult.fPacked[1][2];
		PackedResult.Row2 = Row2Column4;
		Result.Row2W = PackedResult.fPacked[1][3];

		// Third row
		PackedResult.Row3 = Row3Column1;
		Result.Row3X = PackedResult.fPacked[2][0];
		PackedResult.Row3 = Row3Column2;
		Result.Row3Y = PackedResult.fPacked[2][1];
		PackedResult.Row3 = Row3Column3;
		Result.Row3Z = PackedResult.fPacked[2][2];
		PackedResult.Row3 = Row3Column4;
		Result.Row3W = PackedResult.fPacked[2][3];

		PackedResult.Row4 = Row4Column1;
		Result.Row4X = PackedResult.fPacked[3][0];
		PackedResult.Row4 = Row4Column2;
		Result.Row4Y = PackedResult.fPacked[3][1];
		PackedResult.Row4 = Row4Column3;
		Result.Row4Z = PackedResult.fPacked[3][2];
		PackedResult.Row4 = Row4Column4;
		Result.Row4W = PackedResult.fPacked[3][3];

		return Result;
	}

	/*@note: Multiplies Matrix A by Matrix B
	@param A: R-value reference to a 4 by 4 matrix
	@param B: R-value reference to a 4 by 4 matrix
	@return: Resulting product*/
	Matrix4x4 MatrixProduct(Matrix4x4&& A, Matrix4x4&& B)
	{
		Matrix4x4 Result = {};

		// Get the columns
		PackedVector Column1 = LoadVector(B.Row1X, B.Row2X, B.Row3X, B.Row4X);
		PackedVector Column2 = LoadVector(B.Row1Y, B.Row2Y, B.Row3Y, B.Row4Y);
		PackedVector Column3 = LoadVector(B.Row1Z, B.Row2Z, B.Row3Z, B.Row4Z);
		PackedVector Column4 = LoadVector(B.Row1W, B.Row2W, B.Row3W, B.Row4W);

		// And store the dot product of the rows of A with
		// the columns of B to PackedVectors

		/*First row*/
		PackedVector Row1Column1 = DotProduct(A.Row1, Column1, 0xFF);
		PackedVector Row1Column2 = DotProduct(A.Row1, Column2, 0xFF);
		PackedVector Row1Column3 = DotProduct(A.Row1, Column3, 0xFF);
		PackedVector Row1Column4 = DotProduct(A.Row1, Column4, 0xFF);

		/*Second row*/
		PackedVector Row2Column1 = DotProduct(A.Row2, Column1, 0xFF);
		PackedVector Row2Column2 = DotProduct(A.Row2, Column2, 0xFF);
		PackedVector Row2Column3 = DotProduct(A.Row2, Column3, 0xFF);
		PackedVector Row2Column4 = DotProduct(A.Row2, Column4, 0xFF);

		/*Third row*/
		PackedVector Row3Column1 = DotProduct(A.Row3, Column1, 0xFF);
		PackedVector Row3Column2 = DotProduct(A.Row3, Column2, 0xFF);
		PackedVector Row3Column3 = DotProduct(A.Row3, Column3, 0xFF);
		PackedVector Row3Column4 = DotProduct(A.Row3, Column4, 0xFF);

		/*Fourth row*/
		PackedVector Row4Column1 = DotProduct(A.Row4, Column1, 0xFF);
		PackedVector Row4Column2 = DotProduct(A.Row4, Column2, 0xFF);
		PackedVector Row4Column3 = DotProduct(A.Row4, Column3, 0xFF);
		PackedVector Row4Column4 = DotProduct(A.Row4, Column4, 0xFF);

		// Get a packed matrix - union
		SSE::UPackedMatrix PackedResult = {};

		// And store them accordingly - indexing isn't needed but is used to make things clear
		// First row
		PackedResult.Row1 = Row1Column1;
		Result.Row1X = PackedResult.fPacked[0][0];
		PackedResult.Row1 = Row1Column2;
		Result.Row1Y = PackedResult.fPacked[0][1];
		PackedResult.Row1 = Row1Column3;
		Result.Row1Z = PackedResult.fPacked[0][2];
		PackedResult.Row1 = Row1Column4;
		Result.Row1W = PackedResult.fPacked[0][3];

		// Second row
		PackedResult.Row2 = Row2Column1;
		Result.Row2X = PackedResult.fPacked[1][0];
		PackedResult.Row2 = Row2Column2;
		Result.Row2Y = PackedResult.fPacked[1][1];
		PackedResult.Row2 = Row2Column3;
		Result.Row2Z = PackedResult.fPacked[1][2];
		PackedResult.Row2 = Row2Column4;
		Result.Row2W = PackedResult.fPacked[1][3];

		// Third row
		PackedResult.Row3 = Row3Column1;
		Result.Row3X = PackedResult.fPacked[2][0];
		PackedResult.Row3 = Row3Column2;
		Result.Row3Y = PackedResult.fPacked[2][1];
		PackedResult.Row3 = Row3Column3;
		Result.Row3Z = PackedResult.fPacked[2][2];
		PackedResult.Row3 = Row3Column4;
		Result.Row3W = PackedResult.fPacked[2][3];

		PackedResult.Row4 = Row4Column1;
		Result.Row4X = PackedResult.fPacked[3][0];
		PackedResult.Row4 = Row4Column2;
		Result.Row4Y = PackedResult.fPacked[3][1];
		PackedResult.Row4 = Row4Column3;
		Result.Row4Z = PackedResult.fPacked[3][2];
		PackedResult.Row4 = Row4Column4;
		Result.Row4W = PackedResult.fPacked[3][3];

		return Result;
	}

	/*@note: Transforms RowVector by the Matrix Other - 
	@param RowVector: L-value reference to a row vector
	@param Other: L-value reference to a matrix
	@return: Transformed vector*/
	Vector4 TransformRowVector(Vector4& RowVector, Matrix4x4& Other)
	{

		Vector4 RetVal = {};
		// Store the columns
		PackedVector Column1 = LoadVector(Other.Row1X, Other.Row2X, Other.Row3X, Other.Row4X);
		PackedVector Column2 = LoadVector(Other.Row1Y, Other.Row2Y, Other.Row3Y, Other.Row4Y);
		PackedVector Column3 = LoadVector(Other.Row1Z, Other.Row2Z, Other.Row3Z, Other.Row4Z);
		PackedVector Column4 = LoadVector(Other.Row1W, Other.Row2W, Other.Row3W, Other.Row4W);

		// And get the dot product of the rows of the vectors with the columns of the matrix
		PackedVector Row1Column1 = DotProduct(RowVector.xmmVector, Column1, 0xFF);
		PackedVector Row1Column2 = DotProduct(RowVector.xmmVector, Column2, 0xFF);
		PackedVector Row1Column3 = DotProduct(RowVector.xmmVector, Column3, 0xFF);
		PackedVector Row1Column4 = DotProduct(RowVector.xmmVector, Column4, 0xFF);

		// Get a packed vector union
		SSE::UPackedVector PackedUnion = {};

		// And store accordingly - indexing isn't needed but is done for clarity
		PackedUnion.Vector = Row1Column1;
		RetVal.x = PackedUnion.fPacked[0];

		PackedUnion.Vector = Row1Column2;
		RetVal.y = PackedUnion.fPacked[1];

		PackedUnion.Vector = Row1Column3;
		RetVal.z = PackedUnion.fPacked[2];

		PackedUnion.Vector = Row1Column4;
		RetVal.w = PackedUnion.fPacked[3];

		return RetVal;
	}

	/*@note: Transforms RowVector by the Matrix Other -
	@param RowVector: R-value reference to a row vector
	@param Other: L-value reference to a matrix
	@return: Transformed vector*/
	Vector4 TransformRowVector(Vector4&& RowVector, Matrix4x4& Other)
	{

		Vector4 RetVal = {};
		PackedVector Column1 = LoadVector(Other.Row1X, Other.Row2X, Other.Row3X, Other.Row4X);
		PackedVector Column2 = LoadVector(Other.Row1Y, Other.Row2Y, Other.Row3Y, Other.Row4Y);
		PackedVector Column3 = LoadVector(Other.Row1Z, Other.Row2Z, Other.Row3Z, Other.Row4Z);
		PackedVector Column4 = LoadVector(Other.Row1W, Other.Row2W, Other.Row3W, Other.Row4W);

		PackedVector Row1Column1 = DotProduct(RowVector.xmmVector, Column1, 0xFF);
		PackedVector Row1Column2 = DotProduct(RowVector.xmmVector, Column2, 0xFF);
		PackedVector Row1Column3 = DotProduct(RowVector.xmmVector, Column3, 0xFF);
		PackedVector Row1Column4 = DotProduct(RowVector.xmmVector, Column4, 0xFF);

		SSE::UPackedVector PackedUnion = {};

		PackedUnion.Vector = Row1Column1;
		RetVal.x = PackedUnion.fPacked[0];

		PackedUnion.Vector = Row1Column2;
		RetVal.y = PackedUnion.fPacked[1];

		PackedUnion.Vector = Row1Column3;
		RetVal.z = PackedUnion.fPacked[2];

		PackedUnion.Vector = Row1Column4;
		RetVal.w = PackedUnion.fPacked[3];

		return RetVal;
	}
	/*@note: Transforms RowVector by the Matrix Other -
	@param RowVector: L-value reference to a row vector
	@param Other: R-value reference to a matrix
	@return: Transformed vector*/
	Vector4 TransformRowVector(Vector4& RowVector, Matrix4x4&& Other)
	{

		Vector4 RetVal = {};
		PackedVector Column1 = LoadVector(Other.Row1X, Other.Row2X, Other.Row3X, Other.Row4X);
		PackedVector Column2 = LoadVector(Other.Row1Y, Other.Row2Y, Other.Row3Y, Other.Row4Y);
		PackedVector Column3 = LoadVector(Other.Row1Z, Other.Row2Z, Other.Row3Z, Other.Row4Z);
		PackedVector Column4 = LoadVector(Other.Row1W, Other.Row2W, Other.Row3W, Other.Row4W);

		PackedVector Row1Column1 = DotProduct(RowVector.xmmVector, Column1, 0xFF);
		PackedVector Row1Column2 = DotProduct(RowVector.xmmVector, Column2, 0xFF);
		PackedVector Row1Column3 = DotProduct(RowVector.xmmVector, Column3, 0xFF);
		PackedVector Row1Column4 = DotProduct(RowVector.xmmVector, Column4, 0xFF);

		SSE::UPackedVector PackedUnion = {};

		PackedUnion.Vector = Row1Column1;
		RetVal.x = PackedUnion.fPacked[0];

		PackedUnion.Vector = Row1Column2;
		RetVal.y = PackedUnion.fPacked[1];

		PackedUnion.Vector = Row1Column3;
		RetVal.z = PackedUnion.fPacked[2];

		PackedUnion.Vector = Row1Column4;
		RetVal.w = PackedUnion.fPacked[3];

		return RetVal;
	}
	/*@note: Transforms RowVector by the Matrix Other -
	@param RowVector: R-value reference to a row vector
	@param Other: R-value reference to a matrix
	@return: Transformed vector*/
	Vector4 TransformRowVector(Vector4&& RowVector, Matrix4x4&& Other)
	{

		Vector4 RetVal = {};
		PackedVector Column1 = LoadVector(Other.Row1X, Other.Row2X, Other.Row3X, Other.Row4X);
		PackedVector Column2 = LoadVector(Other.Row1Y, Other.Row2Y, Other.Row3Y, Other.Row4Y);
		PackedVector Column3 = LoadVector(Other.Row1Z, Other.Row2Z, Other.Row3Z, Other.Row4Z);
		PackedVector Column4 = LoadVector(Other.Row1W, Other.Row2W, Other.Row3W, Other.Row4W);

		PackedVector Row1Column1 = DotProduct(RowVector.xmmVector, Column1, 0xFF);
		PackedVector Row1Column2 = DotProduct(RowVector.xmmVector, Column2, 0xFF);
		PackedVector Row1Column3 = DotProduct(RowVector.xmmVector, Column3, 0xFF);
		PackedVector Row1Column4 = DotProduct(RowVector.xmmVector, Column4, 0xFF);

		SSE::UPackedVector PackedUnion = {};

		PackedUnion.Vector = Row1Column1;
		RetVal.x = PackedUnion.fPacked[0];

		PackedUnion.Vector = Row1Column2;
		RetVal.y = PackedUnion.fPacked[1];

		PackedUnion.Vector = Row1Column3;
		RetVal.z = PackedUnion.fPacked[2];

		PackedUnion.Vector = Row1Column4;
		RetVal.w = PackedUnion.fPacked[3];

		return RetVal;
	}
	/*@note: Transforms a column vector by this matrix. This is shorter due to the way
	vectors are stored.
	@param Other: L-value reference matrix that transforms ColumnVector
	@param ColumnVector: L-value reference to the column vector that is to be transformed by this matrix
	@return: Transformed vector*/
	Vector4 TransformColumnVector(Matrix4x4& Other, Vector4& ColumnVector)
	{
		Matrix4x4 Temp = {};
		Vector4 RetVal = {};

		// Store the rows in a temporary
		Temp.Row1 = DotProduct(Other.Row1, ColumnVector.xmmVector, 0xFF);
		Temp.Row2 = DotProduct(Other.Row2, ColumnVector.xmmVector, 0xFF);
		Temp.Row3 = DotProduct(Other.Row3, ColumnVector.xmmVector, 0xFF);
		Temp.Row4 = DotProduct(Other.Row4, ColumnVector.xmmVector, 0xFF);
		
		// Get a packed vector union
		SSE::UPackedVector PackedUnion = {};

		// And store accordingly
		PackedUnion.Vector = Temp.Row1;
		RetVal.x = PackedUnion.fPacked[0];

		PackedUnion.Vector = Temp.Row2;
		RetVal.y = PackedUnion.fPacked[1];

		PackedUnion.Vector = Temp.Row3;
		RetVal.z = PackedUnion.fPacked[2];

		PackedUnion.Vector = Temp.Row4;
		RetVal.w = PackedUnion.fPacked[3];

		return RetVal;
	}

	/*@note: Transforms a column vector by this matrix. This is shorter due to the way
	vectors are stored.
	@param Other: R-value reference matrix that transforms ColumnVector
	@param ColumnVector: L-value reference to the column vector that is to be transformed by this matrix
	@return: Transformed vector*/
	Vector4 TransformColumnVector(Matrix4x4&& Other, Vector4& ColumnVector)
	{
		Matrix4x4 Temp = {};
		Vector4 RetVal = {};

		Temp.Row1 = DotProduct(Other.Row1, ColumnVector.xmmVector, 0xFF);
		Temp.Row2 = DotProduct(Other.Row2, ColumnVector.xmmVector, 0xFF);
		Temp.Row3 = DotProduct(Other.Row3, ColumnVector.xmmVector, 0xFF);
		Temp.Row4 = DotProduct(Other.Row4, ColumnVector.xmmVector, 0xFF);

		SSE::UPackedVector PackedUnion = {};
		PackedUnion.Vector = Temp.Row1;
		RetVal.x = PackedUnion.fPacked[0];

		PackedUnion.Vector = Temp.Row2;
		RetVal.y = PackedUnion.fPacked[1];

		PackedUnion.Vector = Temp.Row3;
		RetVal.z = PackedUnion.fPacked[2];

		PackedUnion.Vector = Temp.Row4;
		RetVal.w = PackedUnion.fPacked[3];

		return RetVal;
	}
	/*@note: Transforms a column vector by this matrix. This is shorter due to the way
	vectors are stored.
	@param Other: L-value reference matrix that transforms ColumnVector
	@param ColumnVector: R-value reference to the column vector that is to be transformed by this matrix
	@return: Transformed vector*/
	Vector4 TransformColumnVector(Matrix4x4& Other, Vector4&& ColumnVector)
	{
		Matrix4x4 Temp = {};
		Vector4 RetVal = {};

		Temp.Row1 = DotProduct(Other.Row1, ColumnVector.xmmVector, 0xFF);
		Temp.Row2 = DotProduct(Other.Row2, ColumnVector.xmmVector, 0xFF);
		Temp.Row3 = DotProduct(Other.Row3, ColumnVector.xmmVector, 0xFF);
		Temp.Row4 = DotProduct(Other.Row4, ColumnVector.xmmVector, 0xFF);

		SSE::UPackedVector PackedUnion = {};
		PackedUnion.Vector = Temp.Row1;
		RetVal.x = PackedUnion.fPacked[0];

		PackedUnion.Vector = Temp.Row2;
		RetVal.y = PackedUnion.fPacked[1];

		PackedUnion.Vector = Temp.Row3;
		RetVal.z = PackedUnion.fPacked[2];

		PackedUnion.Vector = Temp.Row4;
		RetVal.w = PackedUnion.fPacked[3];

		return RetVal;
	}
	/*@note: Transforms a column vector by this matrix. This is shorter due to the way
	vectors are stored.
	@param Other: R-value reference matrix that transforms ColumnVector
	@param ColumnVector: R-value reference to the column vector that is to be transformed by this matrix
	@return: Transformed vector*/
	Vector4 TransformColumnVector(Matrix4x4&& Other, Vector4&& ColumnVector)
	{
		Matrix4x4 Temp = {};
		Vector4 RetVal = {};

		Temp.Row1 = DotProduct(Other.Row1, ColumnVector.xmmVector, 0xFF);
		Temp.Row2 = DotProduct(Other.Row2, ColumnVector.xmmVector, 0xFF);
		Temp.Row3 = DotProduct(Other.Row3, ColumnVector.xmmVector, 0xFF);
		Temp.Row4 = DotProduct(Other.Row4, ColumnVector.xmmVector, 0xFF);

		SSE::UPackedVector PackedUnion = {};
		PackedUnion.Vector = Temp.Row1;
		RetVal.x = PackedUnion.fPacked[0];

		PackedUnion.Vector = Temp.Row2;
		RetVal.y = PackedUnion.fPacked[1];

		PackedUnion.Vector = Temp.Row3;
		RetVal.z = PackedUnion.fPacked[2];

		PackedUnion.Vector = Temp.Row4;
		RetVal.w = PackedUnion.fPacked[3];

		return RetVal;
	}
	/*@note: Takes the inverse of A using Cramer's rule
	@param Matrix: L-value reference to the matrix */
	void Inverse(Matrix4x4& Matrix)
	{
		Float4x4& TempMatrix = *((Float4x4*)Matrix.MatrixArray);
		Float4x4 Result;
		Float4x4 Minors;
		float pDeterminant[4];

		Minors[0][0] = TempMatrix[2][2] * TempMatrix[3][3] - TempMatrix[2][3] * TempMatrix[3][2];
		Minors[0][1] = TempMatrix[1][2] * TempMatrix[3][3] - TempMatrix[1][3] * TempMatrix[3][2];
		Minors[0][2] = TempMatrix[1][2] * TempMatrix[2][3] - TempMatrix[1][3] * TempMatrix[2][2];

		Minors[1][0] = TempMatrix[2][2] * TempMatrix[3][3] - TempMatrix[2][3] * TempMatrix[3][2];
		Minors[1][1] = TempMatrix[0][2] * TempMatrix[3][3] - TempMatrix[0][3] * TempMatrix[3][2];
		Minors[1][2] = TempMatrix[0][2] * TempMatrix[2][3] - TempMatrix[0][3] * TempMatrix[2][2];

		Minors[2][0] = TempMatrix[1][2] * TempMatrix[3][3] - TempMatrix[1][3] * TempMatrix[3][2];
		Minors[2][1] = TempMatrix[0][2] * TempMatrix[3][3] - TempMatrix[0][3] * TempMatrix[3][2];
		Minors[2][2] = TempMatrix[0][2] * TempMatrix[1][3] - TempMatrix[0][3] * TempMatrix[1][2];

		Minors[3][0] = TempMatrix[1][2] * TempMatrix[2][3] - TempMatrix[1][3] * TempMatrix[2][2];
		Minors[3][1] = TempMatrix[0][2] * TempMatrix[2][3] - TempMatrix[0][3] * TempMatrix[2][2];
		Minors[3][2] = TempMatrix[0][2] * TempMatrix[1][3] - TempMatrix[0][3] * TempMatrix[1][2];

		pDeterminant[0] = TempMatrix[1][1] * Minors[0][0] - TempMatrix[2][1] * Minors[0][1] + TempMatrix[3][1] * Minors[0][2];
		pDeterminant[1] = TempMatrix[0][1] * Minors[1][0] - TempMatrix[2][1] * Minors[1][1] + TempMatrix[3][1] * Minors[1][2];
		pDeterminant[2] = TempMatrix[0][1] * Minors[2][0] - TempMatrix[1][1] * Minors[2][1] + TempMatrix[3][1] * Minors[2][2];
		pDeterminant[3] = TempMatrix[0][1] * Minors[3][0] - TempMatrix[1][1] * Minors[3][1] + TempMatrix[2][1] * Minors[3][2];

		float Determinant = TempMatrix[0][0] * pDeterminant[0] - TempMatrix[1][0] * pDeterminant[1] + TempMatrix[2][0] * pDeterminant[2] - TempMatrix[3][0] * pDeterminant[3];

		// Ensure matrix is invertible!
		assert(!(MathUtility::IsEqual(Determinant, 0)));
		const float	RDet = 1.0f / Determinant;

		Result[0][0] = RDet * pDeterminant[0];
		Result[0][1] = -RDet * pDeterminant[1];
		Result[0][2] = RDet * pDeterminant[2];
		Result[0][3] = -RDet * pDeterminant[3];
		Result[1][0] = -RDet * (TempMatrix[1][0] * Minors[0][0] - TempMatrix[2][0] * Minors[0][1] + TempMatrix[3][0] * Minors[0][2]);
		Result[1][1] = RDet * (TempMatrix[0][0] * Minors[1][0] - TempMatrix[2][0] * Minors[1][1] + TempMatrix[3][0] * Minors[1][2]);
		Result[1][2] = -RDet * (TempMatrix[0][0] * Minors[2][0] - TempMatrix[1][0] * Minors[2][1] + TempMatrix[3][0] * Minors[2][2]);
		Result[1][3] = RDet * (TempMatrix[0][0] * Minors[3][0] - TempMatrix[1][0] * Minors[3][1] + TempMatrix[2][0] * Minors[3][2]);
		Result[2][0] = RDet * (
			TempMatrix[1][0] * (TempMatrix[2][1] * TempMatrix[3][3] - TempMatrix[2][3] * TempMatrix[3][1]) -
			TempMatrix[2][0] * (TempMatrix[1][1] * TempMatrix[3][3] - TempMatrix[1][3] * TempMatrix[3][1]) +
			TempMatrix[3][0] * (TempMatrix[1][1] * TempMatrix[2][3] - TempMatrix[1][3] * TempMatrix[2][1])
			);
		Result[2][1] = -RDet * (
			TempMatrix[0][0] * (TempMatrix[2][1] * TempMatrix[3][3] - TempMatrix[2][3] * TempMatrix[3][1]) -
			TempMatrix[2][0] * (TempMatrix[0][1] * TempMatrix[3][3] - TempMatrix[0][3] * TempMatrix[3][1]) +
			TempMatrix[3][0] * (TempMatrix[0][1] * TempMatrix[2][3] - TempMatrix[0][3] * TempMatrix[2][1])
			);
		Result[2][2] = RDet * (
			TempMatrix[0][0] * (TempMatrix[1][1] * TempMatrix[3][3] - TempMatrix[1][3] * TempMatrix[3][1]) -
			TempMatrix[1][0] * (TempMatrix[0][1] * TempMatrix[3][3] - TempMatrix[0][3] * TempMatrix[3][1]) +
			TempMatrix[3][0] * (TempMatrix[0][1] * TempMatrix[1][3] - TempMatrix[0][3] * TempMatrix[1][1])
			);
		Result[2][3] = -RDet * (
			TempMatrix[0][0] * (TempMatrix[1][1] * TempMatrix[2][3] - TempMatrix[1][3] * TempMatrix[2][1]) -
			TempMatrix[1][0] * (TempMatrix[0][1] * TempMatrix[2][3] - TempMatrix[0][3] * TempMatrix[2][1]) +
			TempMatrix[2][0] * (TempMatrix[0][1] * TempMatrix[1][3] - TempMatrix[0][3] * TempMatrix[1][1])
			);
		Result[3][0] = -RDet * (
			TempMatrix[1][0] * (TempMatrix[2][1] * TempMatrix[3][2] - TempMatrix[2][2] * TempMatrix[3][1]) -
			TempMatrix[2][0] * (TempMatrix[1][1] * TempMatrix[3][2] - TempMatrix[1][2] * TempMatrix[3][1]) +
			TempMatrix[3][0] * (TempMatrix[1][1] * TempMatrix[2][2] - TempMatrix[1][2] * TempMatrix[2][1])
			);
		Result[3][1] = RDet * (
			TempMatrix[0][0] * (TempMatrix[2][1] * TempMatrix[3][2] - TempMatrix[2][2] * TempMatrix[3][1]) -
			TempMatrix[2][0] * (TempMatrix[0][1] * TempMatrix[3][2] - TempMatrix[0][2] * TempMatrix[3][1]) +
			TempMatrix[3][0] * (TempMatrix[0][1] * TempMatrix[2][2] - TempMatrix[0][2] * TempMatrix[2][1])
			);
		Result[3][2] = -RDet * (
			TempMatrix[0][0] * (TempMatrix[1][1] * TempMatrix[3][2] - TempMatrix[1][2] * TempMatrix[3][1]) -
			TempMatrix[1][0] * (TempMatrix[0][1] * TempMatrix[3][2] - TempMatrix[0][2] * TempMatrix[3][1]) +
			TempMatrix[3][0] * (TempMatrix[0][1] * TempMatrix[1][2] - TempMatrix[0][2] * TempMatrix[1][1])
			);
		Result[3][3] = RDet * (
			TempMatrix[0][0] * (TempMatrix[1][1] * TempMatrix[2][2] - TempMatrix[1][2] * TempMatrix[2][1]) -
			TempMatrix[1][0] * (TempMatrix[0][1] * TempMatrix[2][2] - TempMatrix[0][2] * TempMatrix[2][1]) +
			TempMatrix[2][0] * (TempMatrix[0][1] * TempMatrix[1][2] - TempMatrix[0][2] * TempMatrix[1][1])
			);

		memcpy(Matrix.MatrixArray, &Result, sizeof(Float4x4));
	}
}