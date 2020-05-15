#include <iostream>
#include "isolPnt.h"
#include "poly2.h"

extern "C" {
#include "interval.h"
}

using namespace std;

static vector<int> testCrossings(double* timing, const isolPnt& pt, const poly2& C1, const poly2& C2, bool horiz, bool useMax = false) {
    polynom x; x[1] = 1;
    intvl i = horiz ? pt.y : pt.x;
    intvl I = horiz ? pt.x : pt.y;
    rational v = useMax ? i.vMax() : i.vMin();
    polynom V(v);
    polynom P1 = horiz ? C1(x, V) : C1(V, x);
    polynom P2 = horiz ? C2(x, V) : C2(V, x);
   double timeIsol = get_cpu_time();
    auto rts1 = polySlv(P1.intForm()).roots();
    auto rts2 = polySlv(P2.intForm()).roots();
    auto rtsA = polySlv((P1*P2).intForm()).roots();
   timeIsol = get_cpu_time() - timeIsol;
   *timing += timeIsol;
    int ct1 = rts1.size();
    int ct2 = rts2.size();
    int ctT = rtsA.size();
    long k = I.logDen();
    for(intvl& r : rtsA)
        if(r.logDen() > k) k = r.logDen();
    for(intvl& r : rts1) r.refine(k);
    for(intvl& r : rts2) r.refine(k);
    vector<int> crosses(ctT);
    int i1 = 0, i2 = 0;
    while(i1 < ct1 && rts1[i1].vMax() <= I.vMin()) ++i1;
    while(i2 < ct2 && rts2[i2].vMax() <= I.vMin()) ++i2;
    for(int idx = 0; idx < ctT; ++idx) {
        crosses[idx] = 0;
        intvl& r = rtsA[idx];
        if(r.vMax() <= I.vMin() || r.vMin() >= I.vMax()) continue;
        if(i1 < ct1 && rts1[i1].vMin() >= r.vMin()
                    && rts1[i1].vMax() <= r.vMax()) {
            crosses[idx] += 1; ++i1;
        }
        if(i2 < ct2 && rts2[i2].vMin() >= r.vMin()
                    && rts2[i2].vMax() <= r.vMax()) {
            crosses[idx] += 2; ++i2;
        }
        //TODO: remove non-(v=3) glance crossings
    }
    return crosses;
}

static void pushVec(vector<int>& dst, const vector<int>& src, bool reverse = false) {
    if(reverse)
        dst.insert(dst.end(), src.crbegin(), src.crend());
    else
        dst.insert(dst.end(), src.cbegin(), src.cend());
}

static int testPointFixed(double* timing, const isolPnt& pt, const poly2& C1, const poly2& C2) {
    bool exX = pt.x.isExact();
    bool exY = pt.y.isExact();
    if(exX && exY) {
        if(C1(pt.x.vMin(), pt.y.vMin()) != rational(0)) return -1;
        if(C2(pt.x.vMin(), pt.y.vMin()) != rational(0)) return -1;
        return 1;
    }
    if(exX || exY) {
        vector<int> cr = testCrossings(timing, pt, C1, C2, exY);
        for(int v : cr) if(v == 3) return 1;
        return -1;
    }

    vector<int> cr;
    pushVec(cr, testCrossings(timing, pt, C1, C2, false, false), false);
    pushVec(cr, testCrossings(timing, pt, C1, C2, true , true ), false);
    pushVec(cr, testCrossings(timing, pt, C1, C2, false, true ), true );
    pushVec(cr, testCrossings(timing, pt, C1, C2, true , false), true );

    int c1 = 0, c2 = 0;
    for(int v : cr) {
        switch(v) {
            case 1: ++c1; break;
            case 2: ++c2; break;
            case 3:
                return 1;
        }
    }

    if(c1 == 0 || c2 == 0) return -1;
    if(c1 != 2 || c2 != 2) return 0;

    int p = 0;
    bool alt = true;
    for(int v : cr) if(v != 0) {
        if(v == p) {
            alt = false;
            break;
        } else {
            p = v;
        }
    }
    if(alt) {
        return 1;
    }

    return -1; //TODO
}

static bool testPoint(double* timing, isolPnt& pt, const poly2& C1, const poly2& C2) {
    int res = testPointFixed(timing, pt, C1, C2);
    long k = 1;
    while(res == 0) {
        pt.refine(k);
        k *= 2;
        res = testPointFixed(timing, pt, C1, C2);
    }
    return res > 0;
}

void isolPnt::refine(long minK) {
    x.refine(minK);
    y.refine(minK);
}

vector<isolPnt> cci(const poly2& C1, const poly2& C2, double* timings) {
    double timeRes, timeIsol, timeTotal;
     timeTotal = get_cpu_time();
      polynom resX, resY;
      vector<intvl> rtsX, rtsY;
       timeRes = get_cpu_time();
        resX = resultantX(C1, C2);
        resY = resultantY(C1, C2);
       timeRes = get_cpu_time() - timeRes;
      polySlv psX(resY.intForm());
      polySlv psY(resX.intForm());
       timeIsol = get_cpu_time();
        rtsX = psX.roots();
        rtsY = psY.roots();
       timeIsol = get_cpu_time() - timeIsol;
      vector<isolPnt> ret;
      for(size_t i = 0; i < rtsX.size(); ++i) for(size_t j = 0; j < rtsY.size(); ++j) {
          isolPnt pt;
          pt.x = rtsX[i];
          pt.y = rtsY[j];
          if(testPoint(&timeIsol, pt, C1, C2))
              ret.push_back(pt);
      }
     timeTotal = get_cpu_time() - timeTotal;
    if(timings) {
        timings[0] += timeTotal;
        timings[1] += timeIsol;
        timings[2] += timeRes;
    }
    return ret;
}
