#include "poly2.h"

using namespace std;

poly2::poly2(const polynom& c) : coef(1, c) {}

polynom poly2::operator[](size_t idx) const {
    return idx >= coef.size() ? polynom(0) : coef[idx];
}
polynom& poly2::operator[](size_t idx) {
    if(idx >= coef.size()) coef.resize(idx + 1, polynom(0));
    return coef[idx];
}
void poly2::reduceSize() {
    coef.resize(degreeY() + 1);
    for(polynom& p : coef) p.reduceSize();
}

void poly2::swapXY() {
    size_t degX = degreeX();
    size_t degY = degreeY();
    size_t deg = degX > degY ? degX : degY;
    for(size_t i = 0; i <= deg; ++i)
        for(size_t j = 0; j < i; ++j)
            swap((*this)[i][j], (*this)[j][i]);
    reduceSize();
}

size_t poly2::degreeX() const {
    size_t deg = 0;
    for(const polynom& p : coef) {
        size_t d = p.degree();
        if(d > deg) deg = d;
    }
    return deg;
}

size_t poly2::degreeY() const {
    size_t deg = coef.size() - 1;
    while(deg && (*this)[deg].degree() == 0 && (*this)[deg][0] == 0) --deg;
    return deg;
}

static void printTerm(ostream& out, rational v, size_t nx, size_t ny, bool isFirst) {
    if(v == 0) return;
    if(!isFirst) out << (v > 0 ? " + " : " - ");
    else if(v < 0) out << "-";
    if(v < 0) v = -v;
    if(v != 1 || nx + ny == 0) out << v;
    switch(nx) {
        case 0: break;
        case 1: out << "x"; break;
        default: out << "x^{" << nx << "}"; break;
    }
    switch(ny) {
        case 0: break;
        case 1: out << "y"; break;
        default: out << "y^{" << ny << "}"; break;
    }
}

ostream& operator<<(ostream& out, const poly2& p) {
    size_t degX = p.degreeX();
    size_t degY = p.degreeY();
    bool first = true;
    if(degX == 0 && degY == 0) out << p[0][0];
    else for(size_t i = degX + 1; i;) {
        --i;
        for(size_t j = degY + 1; j;) {
            --j;
            printTerm(out, p[j][i], i, j, first);
            if(first && p[j][i] != 0) first = false;
        }
    }
    return out;
}

polynom resultantY(poly2 A, poly2 B) {
    size_t degA = A.degreeY();
    size_t degB = B.degreeY();

    // Create the Sylvester matrix M
    size_t n = degA + degB;
    std::vector<std::vector<polynom>> M(n, std::vector<polynom>(n));
    for(size_t i = 0; i < degB; ++i)
        for(size_t j = 0; j <= degA; ++j)
            M[i][j + i] = A[j];
    for(size_t i = 0; i < degA; ++i)
        for(size_t j = 0; j <= degB; ++j)
            M[degB + i][j + i] = B[j];

    // Compute the determinant of M via the Bareiss Algorithm
    bool neg = false;
    for(size_t i, j, k = 0; k < n - 1; ++k) {
        for(i = k; i < n; ++i) if(M[i][k].degree() > 0 || M[i][k][0] != 0) {
            neg ^= i != k;
            std::swap(M[i], M[k]);
            break;
        }
        if(i == n) return polynom(0);
        for(i = k + 1; i < n; ++i) for(j = k + 1; j < n; ++j) {
            M[i][j] = M[i][j]*M[k][k] - M[i][k]*M[k][j];
            if(k) M[i][j] /= M[k - 1][k - 1];
        }
    }
    return M[n - 1][n - 1]*(neg ? -1 : 1);
}

polynom resultantX(poly2 A, poly2 B) {
    A.swapXY();
    B.swapXY();
    return resultantY(A, B);
}
