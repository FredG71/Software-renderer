#ifdef __cplusplus
#include <iostream>

#define Breakpoint() __asm int 3
#define assertf( Condition, Message )												\
if (!(Condition)) {                                                                 \
    std::cerr << " Assertion failed " #Condition " failed in " << __FILE__			\
    << " line " << __LINE__ << ": " << Message << std::endl;						\
    Breakpoint()}
#endif
#define IsAligned( Data, Alignment ) assertf( !(((int32_t)Data % Alignment)), "Data is not aligned on " << Alignment );