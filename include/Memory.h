#pragma once
#include <cstdint>
#include <cstdlib>
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include "Debug.h"

template<typename Type> class IGMalloc
{
	class GChunk;
public:
	inline virtual void* Malloc(uint32_t nSize)											= 0;
	inline virtual void Free(void* inPtr)												= 0;
	virtual void* MemMove(void* _Dst, void* _Src, uint32_t nNum)						= 0;
	virtual void* MemCopy(void* __restrict _Dst, void*__restrict _Src, uint32_t nNum)	= 0;
	virtual void* Alloca(size_t nSize)													= 0;
	virtual void* MemSet(void* _Ptr, int32_t nValue, uint32_t nNum)						= 0;
	virtual int32_t MemCompare(const void* _Ptr, const void* _Ptr2, uint32_t nCount)	= 0;
	virtual const void* MemChr(const void* _Ptr, int32_t nValue, uint32_t nNum)			= 0;
	virtual void* MemChr(void* _Ptr, int32_t nValue, uint32_t nNum)						= 0;

protected:
	IGMalloc()
		: nTotalPoolSize(0), nFreePoolSize(0)
	{}

	uint32_t nTotalPoolSize;
	uint32_t nFreePoolSize;
};

template< typename Type > class GMalloc : public IGMalloc<Type>
{
public:
	inline void* Malloc(uint32_t nSize);
	inline void Free(void* inPtr);
	void* MemMove(void* _Dst, void* _Src, uint32_t nNum);
	void* MemCopy(void* __restrict _Dst, void*__restrict _Src, uint32_t nNum);
	void* Alloca(size_t nSize);
	void* MemSet(void* _Ptr, int32_t nValue, uint32_t nNum);
	int32_t MemCompare(const void* _Ptr, const void* _Ptr2, uint32_t nCount);
	const void* MemChr(const void* _Ptr, int32_t nValue, uint32_t nNum);
	void* MemChr(void* _Ptr, int32_t nValue, uint32_t nNum);

	static const int8_t  byMinFreeBlockSize = 10;
	// nSize in bytes
	explicit GMalloc(uint32_t nSize)
	{

		pMemPool = static_cast<int8_t*>(_aligned_malloc(nSize, 0x10));
		IsAligned(pMemPool, 0x10);
		nFreePoolSize = nSize - sizeof(GChunk);
		nTotalPoolSize = nSize;

		GChunk NewChunk = GChunk(nSize - sizeof(GChunk));
		NewChunk.Add(pMemPool);
	}
	GMalloc(uint32_t nSize, uint32_t nAlignment)
	{

		pMemPool = static_cast<int8_t*>(_aligned_malloc(nSize, nAlignment));
		IsAligned(pMemPool, nAlignment);
		nFreePoolSize = nSize - sizeof(GChunk);
		nTotalPoolSize = nSize;

		GChunk NewChunk = GChunk(nSize - sizeof(GChunk));
		NewChunk.Add(pMemPool);
	}
	~GMalloc()
	{
		_aligned_free(pMemPool);
		pMemPool = (int8_t*)0x0;
	}
private:
	class GChunk
	{
	public:
		GChunk()
		{}
		GChunk(uint32_t __nSize)
			:	Next(NULL),
				Prev(NULL),
				nUserSize(__nSize),
				bFree(true)
		{}
		void Add(void* _Dst)
		{
			memcpy(_Dst, this, sizeof(GChunk));
		}
		void Get(void* _Src)
		{
			memcpy(this, _Src, sizeof(GChunk));
		}
		GChunk* Next;
		GChunk* Prev;
		uint32_t nUserSize;
		bool bFree;
	};
	int8_t* pMemPool;
};