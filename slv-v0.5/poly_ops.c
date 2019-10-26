/*=============================================================================
  This file is part of SLV/

  Copyright (C) 2015, Elias Tsigaridas (Elias.Tsigaridas@inria.fr)

  SLV is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  SLV is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.

  If interested by a copy of the GNU Lesser General Public License, write to
  the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
  MA 02110-1301, USA. 
=============================================================================*/


#include <math.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>

#include "poly_ops.h"



/* Print the polynomial */
void slv_poly_pretty_print (FILE* out, mpz_t* F, unsigned long deg)
{
    unsigned long i;
    long s;
    
    if ( mpz_sgn(F[0]) != 0 ) {
        mpz_out_str(out, 10, F[0]);
    }
    
    for (i = 1; i <= deg; i++) {
        s = mpz_sgn(F[i]);
        if ( s == 0 ) continue;
        if ( s == 1 ) {
            fprintf(out, "+");
            if ( mpz_cmp_si(F[i], 1) !=  0 ) {
                mpz_out_str( out, 10, F[i]);
            }
        } else {
            mpz_out_str( out, 10, F[i]);
        }
        if ( mpz_cmp_si(F[i], 1) ==  0 ) { 
            fprintf( out, "");
        } else {
            fprintf( out, "*");
        }
        if (i == 1) {
            fprintf( out, "x");
        } else {
            fprintf( out, "x^%ld", i);
        }
    }
    fprintf( out, "\n");
}



void poly_print_list(FILE* out, mpz_t* F, unsigned long deg)
{
    unsigned long i;
    
    for (i = 0; i <= deg; i++) {
        mpz_out_str(out, 10, F[i]); 
        fprintf(out, " ");
    }
    fprintf(out, "\n");
    return;
    
}


mpz_t* slv_poly_read_from_file(FILE *file_in, const unsigned long rev_flag, unsigned long* deg)
{
	mpz_t *P;
	int i;

	fscanf(file_in, "%lu\n", deg);

	P = (mpz_t*) malloc((*deg + 1) * sizeof(mpz_t));

	if ( rev_flag ) {
		for ( i = *deg; i >= 0; i-- )	{
			mpz_init(P[i]);
			mpz_inp_str(P[i], file_in, 10);
		}
	} else {
		for ( i = 0; i <= *deg; i++ )	{
			mpz_init(P[i]);
			mpz_inp_str(P[i], file_in, 10);
		}
	}
    
	return P;
}


inline 
void slv_poly_reverse (mpz_t* F, unsigned long deg)
{
    long i;
    for ( i = 0; i < (deg+1)/2; i++) {
        mpz_swap(F[i], F[deg-i]);
    }
}



/* Returns the sign of F(1/2) */
inline
int slv_poly_sgn_eval_at_half (mpz_t* F, unsigned long deg)
{
	long j, p;
	int ret;
	mpz_t x, y;

	mpz_init(y);
	mpz_init_set_ui(x, 0);

	p = deg + 1;
	for ( j = 0; j <= deg; j++, p-- ) {
		mpz_mul_2exp(y, F[j], p);
		mpz_add(x, x, y);
	}

	mpz_clear(y);
	ret = mpz_sgn(x);
	mpz_clear(x);
	return ret;
}

inline 
int slv_poly_sgn_eval_at_c(mpz_t* P, long d, mpz_t c)
{
    mpz_t r;
    long j;
    /* fprintf(stderr, "c: ");mpz_out_str(stderr, 10, c); fprintf(stderr, " "); */
    mpz_init_set_ui(r, 0);
    for (j = d; j >= 1; j--) {    
         mpz_add(r, r, P[j]);
         mpz_mul(r, r, c);
    }
    mpz_add(r, r, P[0]);
    /* fprintf(stderr, "r: ");mpz_out_str(stderr, 10, r); fprintf(stderr, " "); */
    
    j = mpz_sgn(r);
    mpz_clear(r);

    return j;
}

inline 
int slv_poly_sgn_eval_at_c_2exp(mpz_t* P, long d, mpz_t c, long k)
{
    mpz_t r;
    mpz_t y;
    long j;
    
    mpz_init_set_ui(r, 1);
    mpz_init_set_ui(y, 1);

    /* fprintf(stderr, "c: ");mpz_out_str(stderr, 10, c); fprintf(stderr, " "); */
    /* debug("k %ld", k); */
    
    if (k <= 0) {
        j = slv_poly_sgn_eval_at_c(P, d, c);
        mpz_clear(r);
        mpz_clear(y);
        return j;
    }

    mpz_mul(r, P[d], c);
    for (j = d-1; j >= 1; j--) {    
        mpz_mul_2exp(y, P[j], (d-j)*k);
        mpz_add(r, r, y);
        mpz_mul(r, r, c);
    }
    mpz_mul_2exp(y, P[0], d*k);
    mpz_add(r, r, y);
    
    j = mpz_sgn(r);
    /* fprintf(stderr, "r: ");mpz_out_str(stderr, 10, r); fprintf(stderr, " "); */
  
    mpz_clear(r);
    mpz_clear(y);
    
    return j;
}



/*! \brief Remove (some) content from a polynomial.
 *
 *
 *   Compute the biggest power of two that is a gcd of all the coefficients 
 *   of the input polynomials and divide the coefficients with this nbber.
 *
 *  \param F polynomial with integer coefficients.
 *  \param deg the degree of the input polynomial.
 *  \return The power of the content.
 */
inline
int slv_poly_remove_content_2exp(mpz_t* F, unsigned long deg)
{
	unsigned long cont, i, z;

	i = 0; while ( mpz_sgn( F[i]) == 0 ) i++;
	cont = get_first_1( F[i]);

	for( ; (i <= deg) && cont; i++ ) {
		if ( mpz_sgn( F[i]) != 0 ) {
			z = get_first_1( F[i]);
			if ( z < cont ) cont = z;
		}
	}
	if ( cont == 0 ) return 0;

	for ( i = 0; i <= deg; i++   )
		mpz_fdiv_q_2exp( F[i], F[i], cont);

	return cont;
}

inline 
long slv_poly_root_upper_bound_2exp(mpz_t* F, unsigned long deg)
{
	long q1, q2, p, i, j;
	unsigned long d = deg;

	long  ad_sgn = mpz_sgn( F[d]);

	q1 = LONG_MIN;
	for (i = 0; i < d; ++i) {
		if ( (mpz_sgn( F[i])  == ad_sgn) || (mpz_sgn( F[i])  == 0) ) continue;

		q2 = LONG_MAX;
		for (j = i+1; j <= d; ++j) {
			if ( mpz_sgn( F[j]) != ad_sgn ) continue;
			p = mpz_bitsize(F[i]) - mpz_bitsize(F[j]) - 1;
			q2 = slv_MIN(q2, p/(j-i) +2);
		}
		q1 = slv_MAX( q1, q2);
	}
	if ( q1 == LONG_MIN )
		return q1 = -1;
	return q1+1;
}

inline 
long slv_poly_root_lower_bound_2exp(mpz_t* F, unsigned long deg)
{
    long q1, q2, p, i, j;
    unsigned long d = deg;

    printf("F[0]: "); mpz_out_str(stdout, 10, F[0]); printf("\n");
    check_debug( mpz_sgn(F[0]) != 0, "lower bound with tcoeff = 0");
    long  a0_sgn = mpz_sgn( F[0]);

    q1 = LONG_MIN;
    for (i = d; i > 0; --i) {
    	if ( (mpz_sgn(F[i])  == a0_sgn) || (mpz_sgn( F[i])  == 0) ) continue;

    	q2 = LONG_MAX;
    	for (j = i-1; j >=0; --j) {
    		if ( mpz_sgn( F[j]) != a0_sgn ) continue;
    		p = mpz_bitsize(F[i]) - mpz_bitsize(F[j]) - 1;
    		q2 = slv_MIN(q2, p/(i-j) +2);
    	}
    	q1 = slv_MAX(q1, q2);
    }
    if ( q1 == LONG_MIN )
    	return q1 = -1;
    return -(q1+1);

    SLV_ERROR exit(-1);
}


inline
void __poly_scale_2exp_pos(mpz_t* F, unsigned long deg, long k)
{
	long i, p;

	p = k;
	for (i = 1; i <= deg; i++, p += k) {
		mpz_mul_2exp( F[i], F[i], p);
	}
	return;
}


inline
void __poly_scale_2exp_neg(mpz_t* F, unsigned long deg, long k)
{
	long i, p;

	p = deg * (-k);
   	for ( i = 0; i < deg ; i++, p += k ) {
		mpz_mul_2exp( F[i], F[i], p);
	}
   	return;
}

inline 
int slv_poly_scale_2exp(mpz_t* F, unsigned long deg, long k)
{
    (k > 0 ? __poly_scale_2exp_pos( F, deg, k) : __poly_scale_2exp_neg( F, deg, k));
	return slv_poly_remove_content_2exp( F, deg);
}



/* Taylor shift by 1. Replaces F by the polynomial F(X+1) */
inline 
void slv_poly_taylor_shift_by_1(mpz_t* F, unsigned long deg)
{
	long i, j;
	for ( i = 0; i <= deg-1; i++ ) {
		for ( j = deg-1 ; j >= i; j-- ) {
			mpz_add( F[j], F[j], F[j+1]);
		}
	}
	return;
}



/* computes the number of sign variations */

long var( mpz_t* f, unsigned long d)
{
	long i, j;
    long v = 0;

    j = 0;
    for ( i=1; i <= d; ++i) {
    	if ( mpz_sgn( f[i]) * mpz_sgn( f[j]) < 0 ) {
    		++v;
    		j = i;
    	}
    }
    return v;
}
