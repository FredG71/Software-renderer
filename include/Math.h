#pragma once
#include <math.h>
#include <cstdio>
#include <cfloat>
#include <smmintrin.h>
#include <intrin.h>
#include <cstdint>
#include "Debug.h"
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
    void Sign(float& fTarget);
    bool IsEqual(float x, float y);
    bool IsEqual(double x, double y);
};

inline bool SSE4Support();

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
    inline Vector4& Normalize();
    inline float Length();
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

    inline Vector3& Normalize();
    inline float Length();
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
        :   Row1(LoadVector(fValue, fValue, fValue, fValue)),
            Row2(LoadVector(fValue, fValue, fValue, fValue)),
            Row3(LoadVector(fValue, fValue, fValue, fValue)),
            Row4(LoadVector(fValue, fValue, fValue, fValue))
    {}
    Matrix4x4()
        :   Row1(LoadVector(1, 0, 0, 0)),
            Row2(LoadVector(0, 1, 0, 0)),
            Row3(LoadVector(0, 0, 1, 0)),
            Row4(LoadVector(0, 0, 0, 1))
    {}
    Matrix4x4(  float fRow1Column1, float fRow1Column2, float fRow1Column3, float fRow1Column4,
                float fRow2Column1, float fRow2Column2, float fRow2Column3, float fRow2Column4,
                float fRow3Column1, float fRow3Column2, float fRow3Column3, float fRow3Column4,
                float fRow4Column1, float fRow4Column2, float fRow4Column3, float fRow4Column4
        )
        :   Row1(LoadVector(fRow1Column1, fRow1Column2, fRow1Column3, fRow1Column4)),
            Row2(LoadVector(fRow2Column1, fRow2Column2, fRow2Column3, fRow2Column4)),
            Row3(LoadVector(fRow3Column1, fRow3Column2, fRow3Column3, fRow3Column4)),
            Row4(LoadVector(fRow4Column1, fRow4Column2, fRow4Column3, fRow4Column4))
    {}

    inline Matrix4x4& operator *(float fScalar);
    __forceinline float Get(uint32_t nRow, uint32_t nColumn);
    Matrix4x4& Inverse(
        EMatrix_Hints Hints,
        EMatrix_Translation_Properties TranslationProperties = EMatrix_Translation_Properties::ERow_Transform,
        EMatrix_Scale_Properties ScaleProperties = EMatrix_Scale_Properties::EUniform_Scale
        );
    inline float GetDeterminant();
    inline Matrix4x4& Transpose();
    Vector4 TransformColumnVector(Vector4&& ColumnVector);
    Vector4 TransformColumnVector(Vector4& ColumnVector);
    Vector4 TransformRowVector(Vector4&& RowVector);
    Vector4 TransformRowVector(Vector4& RowVector);
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
    Matrix4x4 MatrixProduct(Matrix4x4& A, Matrix4x4& B);
    Matrix4x4 MatrixProduct(Matrix4x4&& A, Matrix4x4& B);
    Matrix4x4 MatrixProduct(Matrix4x4& A, Matrix4x4&& B);
    Matrix4x4 MatrixProduct(Matrix4x4&& A, Matrix4x4&& B);
    Vector4 TransformRowVector(Vector4& RowVector, Matrix4x4& Other);
    Vector4 TransformRowVector(Vector4&& RowVector, Matrix4x4& Other);
    Vector4 TransformRowVector(Vector4& RowVector, Matrix4x4&& Other);
    Vector4 TransformRowVector(Vector4&& RowVector, Matrix4x4&& Other);
    Vector4 TransformColumnVector(Matrix4x4& Other, Vector4& ColumnVector);
    Vector4 TransformColumnVector(Matrix4x4&& Other, Vector4& ColumnVector);
    Vector4 TransformColumnVector(Matrix4x4& Other, Vector4&& ColumnVector);
    Vector4 TransformColumnVector(Matrix4x4&& Other, Vector4&& ColumnVector);
    Vector3 ToVector3(Vector4& Target);
    Vector4 ToVector4(Vector3& Target, float w = 0);
    Vector3& ToVector3(Vector3& To, Vector4& From);
    Vector4& ToVector4(Vector4& To, Vector3& From);
    void Inverse(Matrix4x4& Matrix);
}

class Plane{
public:
    Plane(Vector3& Normal, float __D)
        :   NormalVector( Normal ),
            D(__D)
    {}
    Plane(Vector3&& Normal, float __D)
        : NormalVector(Normal),
        D(__D)
    {}
    explicit Plane(Vector4& PlaneEquation)
        :   NormalVector(Math::ToVector3(PlaneEquation)),
            D( PlaneEquation.w )
    {}
    explicit Plane(Vector4&& PlaneEquation)
        : NormalVector(Math::ToVector3(PlaneEquation)),
        D(PlaneEquation.w)
    {}

    Vector3 NormalVector;
    float D;
};