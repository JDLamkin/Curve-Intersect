#include "polynom.h"

using namespace std;

polynom::polynom(const rational& c) : coef(1, c) {}

rational polynom::operator[](size_t idx) const {
    return idx >= coef.size() ? 0 : coef[idx];
}
rational& polynom::operator[](size_t idx) {
    if(idx >= coef.size()) coef.resize(idx + 1, 0);
    return coef[idx];
}
void polynom::reduceSize() {
    coef.resize(degree() + 1);
}

polynom& polynom::operator*=(const rational& o) {
    for(size_t idx = coef.size(); idx;) {
        --idx;
        (*this)[idx] *= o;
    }
    reduceSize();
    return *this;
}
polynom& polynom::operator/=(const rational& o) {
    for(size_t idx = coef.size(); idx;) {
        --idx;
        (*this)[idx] /= o;
    }
    reduceSize();
    return *this;
}
polynom& polynom::operator+=(const polynom& o) {
    for(size_t idx = o.coef.size(); idx;) {
        --idx;
        (*this)[idx] += o[idx];
    }
    reduceSize();
    return *this;
}
polynom& polynom::operator-=(const polynom& o) {
    for(size_t idx = o.coef.size(); idx;) {
        --idx;
        (*this)[idx] -= o[idx];
    }
    reduceSize();
    return *this;
}
polynom& polynom::operator*=(const polynom& o) { return *this = *this * o; }
polynom& polynom::operator/=(const polynom& o) { divide(o); return *this; }
polynom& polynom::operator%=(const polynom& o) { return *this = divide(o); }

polynom polynom::operator-() const {
    polynom r = *this;
    for(rational& val : r.coef) val = -val;
    return r;
}
polynom polynom::operator*(const rational& o) const { polynom p = *this; return p *= o; }
polynom polynom::operator/(const rational& o) const { polynom p = *this; return p /= o; }
polynom polynom::operator+(const polynom& o) const { polynom p = *this; return p += o; }
polynom polynom::operator-(const polynom& o) const { polynom p = *this; return p -= o; }
polynom polynom::operator*(const polynom& o) const {
    polynom r;
    for(size_t idx1 = coef.size(); idx1;) {
        --idx1;
        for(size_t idx2 = o.coef.size(); idx2;) {
            --idx2;
            r[idx1 + idx2] += (*this)[idx1] * o[idx2];
        }
    }
    r.reduceSize();
    return r;
}
polynom polynom::operator/(const polynom& o) const { polynom p = *this; return p /= o; }
polynom polynom::operator%(const polynom& o) const { polynom p = *this; return p.divide(o); }

void polynom::differentiate() {
    if(coef.size() == 0) return;
    for(size_t idx = 1; idx < coef.size(); ++idx)
        coef[idx - 1] = coef[idx] * idx;
    coef.resize(coef.size() - 1);
    reduceSize();
}
polynom polynom::derivative() const { polynom p = *this; p.differentiate(); return p; }

polynom polynom::divide(const polynom& o) {
    size_t oDeg = o.degree();
    rational oLeading = o[oDeg];
    reduceSize();
    if(oDeg == 0) {
        *this /= oLeading;
        return polynom();
    }
    polynom quot;
    // Maintain the value of  *this + o*quot
    while(degree() >= oDeg) {
        size_t pow = degree() - oDeg;
        rational c = coef[degree()] / oLeading;
        quot[pow] += c;
        for(size_t idx = oDeg + 1; idx;) {
            --idx;
            (*this)[idx + pow] -= o[idx] * c;
        }
        reduceSize();
    }
    reduceSize(); quot.reduceSize();
    swap(this->coef, quot.coef);
    return quot; // Now remainder
}

size_t polynom::degree() const {
    size_t deg = coef.size() - 1;
    while(deg && (*this)[deg] == 0) --deg;
    return deg;
}

vector<integer> polynom::intForm() const {
    polynom P = *this;
    P.integize();
    vector<integer> coefI(P.degree() + 1);
    for(size_t i = 0; i < coefI.size(); ++i)
        coefI[i] = P[i].get_num()/P[i].get_den();
    return coefI;
}

void polynom::integize() {
    size_t deg = degree();
    for(size_t i = 0; i <= deg; ++i)
        *this *= (*this)[i].get_den();
}

static void printTerm(ostream& out, rational v, size_t n, bool isFirst) {
    if(v == 0) return;
    if(!isFirst) out << (v > 0 ? " + " : " - ");
    else if(v < 0) out << "-";
    if(v < 0) v = -v;
    if(v != 1 || n == 0) out << v;
    switch(n) {
        case 0: break;
        case 1: out << "x"; break;
        default: out << "x^{" << n << "}"; break;
    }
}

ostream& operator<<(ostream& out, const polynom& p) {
    size_t deg = p.degree();
    if(deg == 0) out << p[0];
    else for(size_t i = deg + 1; i;) {
        --i;
        printTerm(out, p[i], i, i == deg);
    }
    return out;
}
