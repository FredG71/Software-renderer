#pragma once
#include <cstdint>
#define FLAGS_SSE4  0x080000

/*@note: Returns the CPU's features by calling the CPUID instruction.
@return: Feature flags*/
__forceinline uint32_t CPUID_Flags()
{
	uint32_t nFlags = 0;

	__asm{
		push eax
			push ebx
			push ecx
			push edx

			mov eax, 0x1
			cpuid
			mov nFlags, edx

			pop edx
			pop ecx
			pop ebx
			pop eax
	}
	return nFlags;
}

const int32_t nFlags = CPUID_Flags();

/*@return: Whether SSE4 instructions are supported*/
inline bool GetSSE4Support()
{
	return(nFlags & FLAGS_SSE4) != 0;
}