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
/*
 * vca_solver_1.h
 *
 *  Created on: Nov 4, 2013
 *      Author: elias
 */

#ifndef VCA_SOLVER_1_H_
#define VCA_SOLVER_1_H_

#include <gmp.h>

#include "dbg.h"
#include "info.h"
#include "interval.h"
#include "poly_ops.h"

#define TOT_POS 	-1


#define ALL_COEFF_POS   -1

/* Keeps the depth of the tree. */
#define SLV_MAX_DEPTH 2000






void print_root(FILE *stream, slv_bintvl_t z);
void print_roots_all(FILE* stream, slv_bintvl_t* roots, long nbr);



slv_bintvl_t* VCA(mpz_t* A, unsigned long dg, slv_info_ptr info, int algo);
slv_bintvl_t* VCA_0_infinity(mpz_t* A, unsigned long dg, slv_info_ptr info, int algo);
slv_bintvl_t* VCA_all(mpz_t* A, unsigned long dg, slv_info_ptr info_pos, slv_info_ptr info_neg, int algo);

/* Detect if there is one root */
/* Detect the first positive root */
/* Solve in an interval */


/* 
 * Number of sign changes in the coefficients of F(1/(X+1))
 * (Descartes' rule of signs)
 */
int Descartes_test(mpz_t* P, unsigned long deg, long sigh, 
                             long* status,
                             slv_info_ptr info,
                             mpz_t* Q);


void VCA_0_1_mem(mpz_t *FF, 
				unsigned long dg,
                 slv_bintvl_t* vec_bintvl, 
                 slv_info_ptr info);


void slv_adjust_signs(mpz_t* P, long dg, slv_bintvl_t* roots, slv_info_srcptr info);

int slv_refine(mpz_t* P, long dg, slv_bintvl_t I, long t);

/** \brief Refine until the width is <= 2^{-t}
 **/
int slv_refine_until(mpz_t* P, long dg, slv_bintvl_t I, long t);
void slv_refine_until_all(mpz_t* P, long dg, slv_bintvl_t* roots, long nbr, long t);



#endif /* VCA_SOLVER_1_H_ */
