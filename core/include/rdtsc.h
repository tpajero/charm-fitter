#ifndef _RDTSC_H
#define _RDTSC_H

#if defined(__i386__)
static __inline__ unsigned long long rdtsc(void)
{
    unsigned long long int x;
    __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
    return x;
}
#elif defined(__x86_64__)
static __inline__ unsigned long long rdtsc(void)
{
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}
#elif defined(__powerpc__)
static __inline__ unsigned long long rdtsc(void)
{
    unsigned long long int result=0;
    unsigned long int upper, lower,tmp;
    __asm__ volatile(
	    "0:                  \n"
	    "\tmftbu   %0           \n"
	    "\tmftb    %1           \n"
	    "\tmftbu   %2           \n"
	    "\tcmpw    %2,%0        \n"
	    "\tbne     0b         \n"
	    : "=r"(upper),"=r"(lower),"=r"(tmp)
	    );
    result = upper;
    result = result<<32;
    result = result|lower;

    return(result);
}
#elif defined(__aarch64__)
static __inline__ unsigned long long rdtsc(void)
{
    // https://github.com/google/benchmark/blob/v1.1.0/src/cycleclock.h#L116

    // System timer of ARMv8 runs at a different frequency than the CPU's.
    // The frequency is fixed, typically in the range 1-50MHz.  It can be
    // read at CNTFRQ special register.  We assume the OS has set up
    // the virtual timer properly.
    int64_t virtual_timer_value;
    asm volatile("mrs %0, cntvct_el0" : "=r"(virtual_timer_value));
    return virtual_timer_value;
#endif

#endif // _RDTSC_H
