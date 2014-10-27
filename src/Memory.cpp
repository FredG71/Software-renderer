#include "Memory.h"

/*@Note: Allocates a block of memory from the pool
@param nSize: Size in bytes of the memory required
@return: Allocated block*/
template<typename Type> void* GMalloc<Type>::Malloc(uint32_t nSize)
{
	uint32_t nRequiredSize{ nSize + sizeof(GChunk) };
	GChunk* pBlock = (GChunk*)(pMemPool);

	while (pBlock)
	{
		if ( pBlock->bFree && pBlock->nUserSize >= nRequiredSize )
			break;
		pBlock = pBlock->Next;
	}

	int8_t* byBlockData = (int8_t*)pBlock;

	if ( !pBlock )
		return nullptr;

	int32_t nFreeUserSize = pBlock->nUserSize - nRequiredSize;
	if (nFreeUserSize > byMinFreeBlockSize)
	{
		GChunk FreeBlock( nFreeUserSize );
		FreeBlock.Next = pBlock->Next;
		FreeBlock.Prev = pBlock;

		FreeBlock.Add( byBlockData + nRequiredSize );
		if ( FreeBlock.Next )
			FreeBlock.Next->Prev = (GChunk*)( byBlockData + nRequiredSize );
		pBlock->Next = (GChunk*)( byBlockData + nRequiredSize );
		pBlock->nUserSize = nSize;
	}

	nFreePoolSize -= pBlock->nUserSize;
	pBlock->bFree = false;
	return(byBlockData + sizeof(GChunk));
}

/*@Note: Frees a chunk of memory from the pool, ensures that it isn't freed twice
@param _Ptr: Pointer to free*/
template<typename Type> void GMalloc<Type>::Free(void* _Ptr)
{
	if ( !_Ptr )
		return;
	GChunk* pBlock = (GChunk*)((int8_t*)_Ptr - sizeof(GChunk));
	assertf( !pBlock->bFree, "Block is free!");

	if ( pBlock->bFree )
		return;

	int32_t nFullBlockSize = pBlock->nUserSize + sizeof(GChunk);
	nFreePoolSize += pBlock->nUserSize;

	GChunk* Head = pBlock;
	GChunk* pPrev = pBlock->Prev;
	GChunk* pNext = pBlock->Next;

	if (pBlock->Prev && pBlock->Prev->bFree)
	{
		Head = pBlock->Prev;
		pPrev = pBlock->Prev->Prev;
		pNext = pBlock->Next;

		nFullBlockSize += pBlock->Prev->nUserSize + sizeof(GChunk);

		if (pBlock->Next)
		{
			pBlock->Next->Prev = Head;
			if (pBlock->Next->bFree)
			{
				pNext = pBlock->Next->Next;
				if ( pBlock->Next->Next )
					pBlock->Next->Next->Prev = Head;

				nFullBlockSize += pBlock->Next->nUserSize + sizeof(GChunk);
			}
		}
	}
	else
	{
		if (pBlock->Next && pBlock->Next->bFree)
		{
			Head = pBlock;
			pPrev = pBlock->Prev;
			pNext = pBlock->Next->Next;

			nFullBlockSize += pBlock->Next->nUserSize + sizeof(GChunk);
		}
	}
		int8_t* pFreeBlockStart = ( int8_t* )Head;
		uint32_t nFreeUserSize = nFullBlockSize - sizeof(GChunk);

		GChunk FreeBlock( nFreeUserSize );
		FreeBlock.Prev = pPrev;
		FreeBlock.Next = pNext;
		FreeBlock.Add(pFreeBlockStart);
}
/*@Note: Wrapper to memmove
@param _Dst: Destination pointer
@param _Src: Source pointer
@param nNum: Count*/
template<typename Type> void* GMalloc<Type>::MemMove(void* _Dst, void* _Src, uint32_t nNum)
{
	assertf(_Dst && _Src, "_Dst or _Src points to 0x0");
	return memmove( _Dst, _Src, nNum );
}
/*@Note: Wrapper to memcpy
@param _Dst: Destination pointer
@param _Src: Source pointer
@param nNum: Count*/
template<typename Type> void* GMalloc<Type>::MemCopy(void* _Dst, void* _Src, uint32_t nNum)
{
	assertf(_Dst && _Src, "_Dst or _Src points to 0x0");
	return memcpy(_Dst, _Src, nNum);
}
/*@Note: Wrapper to _alloca, allocates memory from the stack, and ensures that nSize is between 0 and 1024,
exceptions are handled
@param nSize: How much to allocate
@return: Allocated memory*/
template<typename Type> void* GMalloc<Type>::Alloca(uint32_t nSize)
{
	int32_t nErrorCode = 0;
	static void* pData = NULL;
	__try
	{
		if (nSize > 0 && nSize < 1024)
			pData = _alloca(nSize);
	}
	__except (GetExceptionCode() == STATUS_STACK_OVERFLOW )
	{
			nErrorCode = _resetstkoflw();
			assertf(!nErrorCode, "Stack couldn't be reset.");
	};
	return pData;
}
/*@Note: Wrapper to memcmp
@param _Ptr: First buffer
@param _Ptr2: Second buffer
@param nCount: Count*/
template<typename Type> int32_t GMalloc<Type>::MemCompare(
	const void* _Ptr,
	const void* _Ptr2,
	uint32_t nCount
	)
{
	assertf(_Ptr && _Ptr2, "_Ptr or _Ptr points to 0x0");
	return memcmp( _Ptr, _Ptr2, nCount );
}
/*@Note: Wrapper to memset
@param _Ptr: Pointer to be set
@param nValue: Value to set
@param nNum: Count*/
template< typename Type > void* GMalloc<Type>::MemSet(void* _Ptr, int32_t nValue, uint32_t nNum)
{
	assertf(_Ptr, "_Ptr points to 0x0.");
	return memset( _Ptr, nValue, nNum );
}
/*@Note: Wrapper to memchr
@param _Ptr: Buffer
@param nValue: Value to look for
@param nNum: Count*/
template<typename Type> const void* GMalloc<Type>::MemChr(const void* _Ptr, int32_t nValue, uint32_t nNum)
{
	assertf(_Ptr, "_Ptr points to 0x0.");
	return memchr( _Ptr, nValue, nNum );
}
/*@Note: Wrapper to memchr
@param _Ptr: Buffer
@param nValue: Value to look for
@param nNum: Count*/
template<typename Type> void* GMalloc<Type>::MemChr(void* _Ptr, int32_t nValue, uint32_t nNum)
{
	assertf(_Ptr, "_Ptr points to 0x0.");
	return memchr(_Ptr, nValue, nNum);
}