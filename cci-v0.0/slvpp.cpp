#include "slvpp.h"

extern "C" {
    #include "dbg.h"

    #include "info.h"
    #include "interval.h"

    #include "poly_ops.h"
    #include "vca_solver_1.h"
};

using namespace std;

static mpq_class getDyadic(mpz_class num, long logDen) {
    if(logDen < 0) logDen = 0;
    mpq_class val;
    val.get_num() = num;
    mpz_ui_pow_ui(val.get_den_mpz_t(), 2, logDen);
    val.canonicalize();
    return val;
}

void polySlv::copy(const polySlv& o) {
    deg = o.deg;
    arr = (mpz_t*)malloc((deg + 1) * sizeof(mpz_t));
    for(unsigned int i = 0; i <= deg; ++i)
        mpz_init_set(arr[i], o.arr[i]);
}
void polySlv::move(polySlv& o) {
    arr = o.arr;
    deg = o.deg;
    o.arr = nullptr;
}
void polySlv::destruct() {
    for(unsigned long i = 0; i <= deg; ++i)
        mpz_clear(arr[i]);
    free(arr);
}
polySlv::polySlv() : polySlv(vector<mpz_class>(1, mpz_class(0))) {}
polySlv::polySlv(const std::vector<mpz_class>& coeffs) {
    deg = coeffs.size() - 1;
    if(coeffs.size() == 0) deg = 0;
    arr = (mpz_t*)malloc((deg + 1) * sizeof(mpz_t));
    for(unsigned int i = 0; i < coeffs.size(); ++i)
        mpz_init_set(arr[i], coeffs[i].get_mpz_t());
    if(coeffs.size() == 0) mpz_init_set_ui(arr[0], 0);
}
vector<intvl> polySlv::roots() const {
    slv_info_t info_pos, info_neg;
    slv_info_init(info_pos);
    slv_info_init(info_neg);

    slv_bintvl_t* rts = VCA_all(arr, deg, info_pos, info_neg, 1);
    vector<intvl> ret(info_pos->nb_roots + info_neg->nb_roots);

    auto F = make_shared<polySlv>(*this);
    for(unsigned long i = 0; i < ret.size(); ++i) {
        ret[i].destruct();
        ret[i].F = F;
        ret[i].copy(rts[i]);
        slv_bintvl_clear(rts[i]);
        ret[i].refine(0);
    }
    free(rts);

    return ret;
}
void intvl::copy(const slv_bintvl* i) {
    interval = (slv_bintvl*)malloc(sizeof(slv_bintvl));
    slv_bintvl_init_set(interval, i);
}
void intvl::copy(const intvl& o) {
    copy(o.interval);
    F = o.F;
}
void intvl::move(intvl& o) {
    interval = o.interval;
    o.interval = nullptr;
    F = o.F;
}
void intvl::destruct() {
    if(interval == nullptr) return;
    slv_bintvl_clear(interval);
    free(interval);
    interval = nullptr;
}
intvl::intvl() : interval(nullptr), F() {}
void intvl::refine(long minK) {
    if(interval == nullptr) return;
	slv_refine_until(F->arr, F->deg, interval, minK);
}
mpq_class intvl::vMin() const {
    return getDyadic(numMin(), logDen());
}
mpq_class intvl::vMax() const {
    return getDyadic(numMax(), logDen());
}
mpq_class intvl::vAve() const {
    return getDyadic(numMin() + numMax(), logDen() + 1);
}
mpz_class intvl::numMin() const {
    if(interval == nullptr) return mpz_class(0);
    return mpz_class(interval->c);//TODO
}
mpz_class intvl::numMax() const {
    return numMin() + ( isExact() ? 0 : 1 );
}
bool intvl::isExact() const {
    if(interval == nullptr) return true;
    return false;//TODO
}
long intvl::logDen() const {
    if(interval == nullptr) return 0;
    return interval->k;//TODO
}
bool intvl::operator<(const intvl& o) const {
    return vMax() <= o.vMin();
}
