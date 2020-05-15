#include <iostream>
#include <iomanip>
#include "isolPnt.h"
#include "poly2.h"

extern "C" {
#include "interval.h"
}

using namespace std;

poly2 getPoly2() {
    poly2 P;
    size_t nx, ny;
    cin >> nx >> ny;
    for(size_t i = 0; i <= ny; ++i)
        for(size_t j = 0; j <= nx; ++j)
            cin >> P[i][j];
    return P;
}

int main0() {
    poly2 A = getPoly2(), B = getPoly2();
    cout << A << " = 0" << endl;
    cout << B << " = 0" << endl;

    double timings[3];
    timings[0] = timings[1] = timings[2] = 0;
    vector<isolPnt> pts = cci(A, B, timings);

    for(isolPnt pt : pts) {
        cout << "max(2^{" << pt.x.logDen() << "} abs(x - " << pt.x.vAve() << "), 2^{" << pt.y.logDen() << "} abs(y - " << pt.y.vAve() << ")) <= 0.5" << endl;
        pt.refine(10);
        cout << "(" << pt.x.vAve() << ", " << pt.y.vAve() << ")" << endl;
    }
    cerr.precision(17);
    cerr << fixed << timings[2] << " " << fixed << timings[1] << endl;

    return 0;
}

void print(const isolPnt& pt) {
    rational x = (pt.x.vMin() + pt.x.vMax())/2;
    rational y = (pt.y.vMin() + pt.y.vMax())/2;
    cout << "(" << mpf_class(x) << ", " << mpf_class(y) << ")" << endl;
}

integer randInt(int L) {
    static gmp_randstate_t state;
    static bool initialized = false;
    if(!initialized) {
        initialized = true;
        gmp_randinit_default(state);
        gmp_randseed_ui(state, (long)(get_cpu_time()*1000000));
    }
    mpz_class neg;
    mpz_class ret;
    mpz_urandomb(neg.get_mpz_t(), state, 1);
    mpz_urandomb(ret.get_mpz_t(), state, L);
    return neg == 0 ? ret : -ret;
}

poly2 makePoly2(int n, int m, int L) {
    poly2 P;
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < m; ++j)
            P[j][i] = randInt(L);
    return P;
}

int main1() {
    int n, L;
    cin >> n >> L;
    ++n;
    int N = 10;
    double timings[3];
    timings[0] = timings[1] = timings[2] = 0;
    for(int idx = 0; idx < N; ++idx) {
        poly2 A = makePoly2(n, n, L), B = makePoly2(n, n, L);
        cci(A, B, timings);
    }
    timings[0] /= N;
    timings[1] /= N;
    timings[2] /= N;
    cerr.precision(9);
    cerr << setw(3) << (n-1) << " "
         << setw(7) << L << " "
         << fixed << timings[0] << " "
         << fixed << timings[1] << " "
         << fixed << timings[2] << " ";
    cerr.precision(8);
    cerr << fixed << setw(11) << 100*timings[1]/timings[0] << " "
         << fixed << setw(11) << 100*timings[2]/timings[0] << endl;
    return 0;
}

int main() {
    return main0();
}
