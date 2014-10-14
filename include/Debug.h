#include <cstdint>

struct RecursionCounter{
	RecursionCounter(uint8_t& nRecurse)
	: __nRecurse__(nRecurse)
		{	__nRecurse__ += 1; }
	~RecursionCounter()
		{	__nRecurse__ -= 1; }
private:
	uint8_t __nRecurse__;
};

#define Recursion()		static uint8_t nCounter;							\
						assert(nCounter == 0);								\
						const RecursionCounter _RecursionCounter(nCounter)

#define ErrorOnEntry() assert(0)