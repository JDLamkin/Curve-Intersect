#ifndef __CPP_UTILS_H__ 
#define __CPP_UTILS_H__ 

#define RAII5(TYPE)                                 \
    TYPE(const TYPE& o) { copy(o); }                \
    TYPE(TYPE&& o) { move(o); }                     \
    ~TYPE() { destruct(); }                         \
    TYPE& operator=(const TYPE& o) {                \
        if(&o != this) {                            \
            destruct();                             \
            copy(o);                                \
        }                                           \
        return *this;                               \
    }                                               \
    TYPE& operator=(TYPE&& o) {                     \
        if(&o != this) {                            \
            destruct();                             \
            move(o);                                \
        }                                           \
        return *this;                               \
    }

#endif
