#ifndef __SLVPP_H__
#define __SLVPP_H__

#include <gmpxx.h>
#include <vector>
#include <memory>
#include "cppUtils.h"

class slv_bintvl;
class intvl;

class polySlv {
private:
    friend class intvl;

    mpz_t* arr;
    unsigned int deg;

    void copy(const polySlv&);
    void move(polySlv&);
    void destruct();

public:
    polySlv();
    explicit polySlv(const std::vector<mpz_class>&);
    RAII5(polySlv)

    std::vector<intvl> roots() const;
};

class intvl {
private:
    friend class polySlv;

    slv_bintvl* interval;
    std::shared_ptr<polySlv> F;

    void copy(const slv_bintvl*);
    void copy(const intvl&);
    void move(intvl&);
    void destruct();

public:
    intvl();
    RAII5(intvl)

    void refine(long minK);

    mpq_class vMin() const;
    mpq_class vMax() const;
    mpq_class vAve() const;
    mpz_class numMin() const;
    mpz_class numMax() const;
    bool isExact() const;
    long logDen() const;

    bool operator<(const intvl&) const;
};
#endif
