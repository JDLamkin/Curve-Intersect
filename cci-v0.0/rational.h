#ifndef _LIBALG_RATIONAL_H_
#define _LIBALG_RATIONAL_H_

#include <gmpxx.h>

typedef mpq_class rational;
typedef mpz_class integer;

/*
static void trigger_div_by_zero() {
    volatile int num = 1;
    volatile int den = 0;
    volatile int qot = num / den;
}*/

#endif
