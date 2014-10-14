#include "Math.h"

template<typename T> T Math::Abs(T Type)
{
	ErrorOnEntry();
	return (T)0x0;
}

template<> float Math::Abs(float fNumber)
{
	__asm{
			finit
			wait
			fld dword ptr[fNumber]
			fabs
			fstp dword ptr[fNumber]
	}
	return fNumber;
}

template<> int32_t Math::Abs(int32_t nNumber)
{
	__asm{
			mov eax, nNumber
			cdq
			xor eax, edx
			sub eax, edx
	}
}