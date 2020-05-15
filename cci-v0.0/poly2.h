#ifndef _LIBALG_POLY2_H_
#define _LIBALG_POLY2_H_

#include "polynom.h"

class poly2 {
private:
    std::vector<polynom> coef;

public:
    poly2(const rational& c = 0) : poly2(polynom(c)) {}
    poly2(const polynom& c);

    polynom operator[](size_t) const;
    polynom& operator[](size_t);
    void reduceSize();

    template<class T>
    T operator()(const T&, const T&) const;

    size_t degreeX() const;
    size_t degreeY() const;
    void swapXY();
};

std::ostream& operator<<(std::ostream&, const poly2&);
polynom resultantX(poly2 A, poly2 B);
polynom resultantY(poly2 A, poly2 B);

template<class T>
T poly2::operator()(const T& x, const T& y) const {
    T sum = T(0);
    T pow = T(1);
    for(size_t i = 0; i < coef.size(); ++i) {
        sum += pow * coef[i](x);
        pow *= y;
    }
    return sum;
}

#endif
