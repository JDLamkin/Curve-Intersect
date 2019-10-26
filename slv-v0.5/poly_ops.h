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


#ifndef __POLY_OPS_H__
#define __POLY_OPS_H__


#include <stdio.h>
#include <stdlib.h>

#include "dbg.h"
#include "interval.h"

#define get_first_1(x) 	mpz_scan1((x), 0)
#define mpz_bitsize(a) 	mpz_sizeinbase((a),2)


#define slv_MIN(a, b) ( ((a) < (b)) ? (a) : (b) )
#define slv_MAX(a, b) ( ((a) > (b)) ? (a) : (b) )


#define SLV_POLY_INIT(R, deg, j)                        \
	R = (mpz_t *) malloc ((deg + 1) * sizeof(mpz_t));   \
	for (j = 0; j <= deg; j++) {                        \
		mpz_init(R[j]);                                 \
    }                                                   \

#define SLV_POLY_INIT_SET(R, F, deg, j)                 \
	R = (mpz_t *) malloc ((deg + 1) * sizeof(mpz_t));   \
	for (j = 0; j <= deg; j++) {                        \
		mpz_init_set(R[j], F[j]);                       \
    }                                                   \

#define SLV_POLY_SET(R, F, deg, j)                      \
	for (j = 0; j <= deg; j++) {                        \
		mpz_set(R[j], F[j]);                            \
    }                                                   \

#define SLV_POLY_CLEAR(R, deg, j)                 \
    for (j = 0; j <= deg; j++) {                  \
		mpz_clear(R[j]);                          \
    }                                             \
    free(R);                                      \


void slv_poly_print_list( FILE* out, mpz_t* F, unsigned long deg);
void slv_poly_pretty_print (FILE* out, mpz_t *F, unsigned long deg);

mpz_t* slv_poly_read_from_file(FILE *file_in, const unsigned long rev_flag, unsigned long* deg);



void slv_poly_reverse (mpz_t* F, unsigned long deg);


int slv_poly_sgn_eval_at_half (mpz_t* F, unsigned long deg);

int slv_poly_sgn_eval_at_c(mpz_t* P, long d, mpz_t c);

int slv_poly_sgn_eval_at_c_2exp(mpz_t* P, long d, mpz_t c, long k);



int slv_poly_remove_content_2exp (mpz_t* F, unsigned long deg);
 

long slv_poly_root_upper_bound_2exp (mpz_t* F, unsigned long deg);

long slv_poly_root_lower_bound_2exp (mpz_t* F, unsigned long deg);

int slv_poly_scale_2exp (mpz_t* F, unsigned long deg, long k);


void slv_poly_taylor_shift_by_1 (mpz_t* F, unsigned long deg);

long var( mpz_t* f, unsigned long d);


#endif /* __POLY_OPS_H__ */

