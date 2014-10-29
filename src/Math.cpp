#include "Math.h"
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

void MathUtility::Sign(float& fTarget)
{
    (*(int32_t*)&fTarget) ^= 0x80000000;
}

bool MathUtility::IsEqual(float x, float y)
{
    return abs(x - y) <= FLT_EPSILON * abs(x);
}

bool MathUtility::IsEqual(double x, double y)
{
    return abs(x - y) <= DBL_EPSILON * abs(x);
}


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


// Normalizes this vector
inline Vector4& Vector4::Normalize()
{
    // Get the square root of the dot product - 4 component
    PackedVector mDot = SquareRoot(DotProduct(xmmVector, xmmVector, 0xFF));
    // divide returned result
    xmmVector = Divide(xmmVector, mDot);
    return *this;
}

// returns the magnitude (length) of this vector
inline float Vector4::Length()
{
    ALIGN_16 float fRetVal[4];
    // Gets the square root of the dot product
    PackedVector mDot = SquareRoot(DotProduct(xmmVector, xmmVector, 0xFF));
    // Store it and return it
    Store(fRetVal, mDot);
    return fRetVal[0];
}


inline Vector3& Vector3::Normalize()
{
    PackedVector mDot = SquareRoot(DotProduct(xmmVector, xmmVector, 0x7F));
    xmmVector = Divide(xmmVector, mDot);
    return *this;
}

inline float Vector3::Length()
{
    ALIGN_16 float fRetVal[4];
    PackedVector mDot = SquareRoot(DotProduct(xmmVector, xmmVector, 0x7F));
    Store(fRetVal, mDot);
    return fRetVal[0];
}



inline Matrix4x4& Matrix4x4::operator *(float fScalar)
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

float Matrix4x4::Get(uint32_t nRow, uint32_t nColumn)
{
    assertf(nRow <= 3 && nRow >= 0 && nColumn <= 3 && nColumn >= 0, "Out of bounds");
    return MatrixArray[nRow][nColumn];
}

/* @note: Computes the inverse of this matrix and gathers a few hints from the caller
in cases where the inversion can be faster, for instance orthogonal matrices etc
@param Hints: Hints to speed up the inversion
@param TranslationProperties: What kind of translation is occuring, is it transforming
a row or column vector?
@param ScaleProperties: Is it a uniform or non-uniform scale?
@return: Inversed matrix*/

Matrix4x4& Matrix4x4::Inverse(
    EMatrix_Hints Hints,
    EMatrix_Translation_Properties TranslationProperties,
    EMatrix_Scale_Properties ScaleProperties
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
                                                            register float fInverse = 1 / MatrixArray[0][0];
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
        assertf(!(MathUtility::IsEqual(Determinant, 0)), "Matrix has no inverse!");
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

inline float Matrix4x4::GetDeterminant()
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
inline Matrix4x4& Matrix4x4::Transpose()
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
Vector4 Matrix4x4::TransformColumnVector(Vector4&& ColumnVector)
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
Vector4 Matrix4x4::TransformColumnVector(Vector4& ColumnVector)
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
Vector4 Matrix4x4::TransformRowVector(Vector4&& RowVector)
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
Vector4 Matrix4x4::TransformRowVector(Vector4& RowVector)
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


/*@note: Multiplies Matrix A by Matrix B
@param A: L-value reference to a 4 by 4 matrix
@param B: L-value reference to a 4 by 4 matrix
@return: Resulting product*/
Matrix4x4 Math::MatrixProduct(Matrix4x4& A, Matrix4x4& B)
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
Matrix4x4 Math::MatrixProduct(Matrix4x4&& A, Matrix4x4& B)
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
Matrix4x4 Math::MatrixProduct(Matrix4x4& A, Matrix4x4&& B)
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
Matrix4x4 Math::MatrixProduct(Matrix4x4&& A, Matrix4x4&& B)
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
Vector4 Math::TransformRowVector(Vector4& RowVector, Matrix4x4& Other)
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
Vector4 Math::TransformRowVector(Vector4&& RowVector, Matrix4x4& Other)
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
Vector4 Math::TransformRowVector(Vector4& RowVector, Matrix4x4&& Other)
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
Vector4 Math::TransformRowVector(Vector4&& RowVector, Matrix4x4&& Other)
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
Vector4 Math::TransformColumnVector(Matrix4x4& Other, Vector4& ColumnVector)
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
Vector4 Math::TransformColumnVector(Matrix4x4&& Other, Vector4& ColumnVector)
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
Vector4 Math::TransformColumnVector(Matrix4x4& Other, Vector4&& ColumnVector)
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
Vector4 Math::TransformColumnVector(Matrix4x4&& Other, Vector4&& ColumnVector)
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
void Math::Inverse(Matrix4x4& Matrix)
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
    assertf(!(MathUtility::IsEqual(Determinant, 0)), "Matrix has no inverse!");
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

Vector3 Math::ToVector3(Vector4& Target)
{
    Vector3 RetVal = {};
    RetVal = *(Vector3*)&Target;
    return RetVal;
}

Vector4 Math::ToVector4(Vector3& Target, float w)
{
    Vector4 RetVal = {};
    RetVal = *(Vector4*)&Target;
    RetVal.w = w;
    return RetVal;
}

Vector3& Math::ToVector3(Vector3& To, Vector4& From)
{
    To = *(Vector3*)&From;
    return To;
}

Vector4& Math::ToVector4(Vector4& To, Vector3& From )
{
    To = *(Vector4*)&From;
    return To;
}

/*
int main()
{
assertf(SSE4Support(), "SSE4 aren't supported on this platform");
Matrix4x4 M = Matrix4x4(1, 2, 3, 4,
5, 6, 7, 8,
9, 10, 11, 12,
13, 14, 15, 16);
Vector4 V = Vector4(2, 32, 4, 89);
V = Math::TransformColumnVector(M, V);
printf("Determinent: %f\n", M.GetDeterminant());

printf("%f, %f, %f, %f\n", V.x, V.y, V.z, V.w);
printf("\n\n\n\n");
printf("%f, %f, %f, %f\n", M.Row1X, M.Row1Y, M.Row1Z, M.Row1W);
printf("%f, %f, %f, %f\n", M.Row2X, M.Row2Y, M.Row2Z, M.Row2W);
printf("%f, %f, %f, %f\n", M.Row3X, M.Row3Y, M.Row3Z, M.Row3W);
printf("%f, %f, %f, %f\n", M.Row4X, M.Row4Y, M.Row4Z, M.Row4W);

M.Inverse(EMatrix_Hints::EUnknown_Matrix, EMatrix_Translation_Properties::ENo_Translation, EMatrix_Scale_Properties::ENo_Scale);

printf("\n\n\n");
printf("%f, %f, %f, %f\n", M.Row1X, M.Row1Y, M.Row1Z, M.Row1W);
printf("%f, %f, %f, %f\n", M.Row2X, M.Row2Y, M.Row2Z, M.Row2W);
printf("%f, %f, %f, %f\n", M.Row3X, M.Row3Y, M.Row3Z, M.Row3W);
printf("%f, %f, %f, %f\n", M.Row4X, M.Row4Y, M.Row4Z, M.Row4W);

getchar();
}*/