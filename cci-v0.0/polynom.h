#ifndef _LIBALG_POLYNOM_H_
#define _LIBALG_POLYNOM_H_

#include <iostream>
#include <vector>
#include "rational.h"

class polynom {
private:
    std::vector<rational> coef;

public:
    polynom(const rational& c = 0);

    rational operator[](size_t) const;
    rational& operator[](size_t);
    void reduceSize();

    polynom& operator*=(const rational&);
    polynom& operator/=(const rational&);
    polynom& operator+=(const polynom&);
    polynom& operator-=(const polynom&);
    polynom& operator*=(const polynom&);
    polynom& operator/=(const polynom&);
    polynom& operator%=(const polynom&);

    polynom operator-() const;
    polynom operator*(const rational&) const;
    polynom operator/(const rational&) const;
    polynom operator+(const polynom&) const;
    polynom operator-(const polynom&) const;
    polynom operator*(const polynom&) const;
    polynom operator/(const polynom&) const;
    polynom operator%(const polynom&) const;

    template<class T>
    T operator()(const T&) const;

    void differentiate();
    polynom derivative() const;

    polynom divide(const polynom&);
    size_t degree() const;

    std::vector<integer> intForm() const;
    void integize();
};

std::ostream& operator<<(std::ostream&, const polynom&);

template<class T>
T polynom::operator()(const T& arg) const {
    T sum = T(0);
    T pow = T(1);
    for(size_t i = 0; i < coef.size(); ++i) {
        sum += pow * coef[i];
        pow *= arg;
    }
    return sum;
}

#endif
