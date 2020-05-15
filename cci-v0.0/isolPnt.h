#ifndef __ISOL_PNT_H__
#define __ISOL_PNT_H__

#include "slvpp.h"

class poly2;

class isolPnt {
public:
    intvl x, y;

    void refine(long minK);
};

std::vector<isolPnt> cci(const poly2& C1, const poly2& C2, double* timings = nullptr);

#endif
