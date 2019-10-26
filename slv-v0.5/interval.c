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
 * interval.c
 *
 *  Created on: Nov 4, 2013
 *      Author: elias
 */

#include <stdlib.h>
#include <time.h>
#include "interval.h"


void slv_bintvl_init( slv_bintvl_ptr a)
{
	a->is_exact = 0;
	a->k 		= 0;
	a->sgn_left = 0;
	a->sign 	= 0;
    
    mpz_init(a->c);

    
}

void slv_bintvl_init_set( slv_bintvl_ptr a, slv_bintvl_srcptr b)
{

	a->is_exact = b->is_exact;
	a->k 		= b->k;
	a->sgn_left = b->sgn_left;
	a->sign 	= b->sign;

	mpz_init_set(a->c, b->c); 
}

void slv_bintvl_clear( slv_bintvl_ptr a)
{
    mpz_clear( a->c);
}


/*! \fn void slv_bintvl_print    (slv_bintvl_srcptr I)
  \brief Prints an interval.
  \param I The interval to print 
*/
void slv_bintvl_print( slv_bintvl_srcptr I)
{

    printf( "I: ");
	mpz_out_str( stdout, 10, I->c);
	printf( "/2^%ld \n", I->k);
}   


void slv_bintvl_get_mid(mpq_t m, slv_bintvl_srcptr I)
{
    mpz_t t;
    mpz_init_set(t, I->c);

    mpq_init(m);

    if ( I->k >= 0 ) {
        mpz_mul_2exp(t, t, 1);
        mpz_add_ui(t, t, 1);
      
        mpq_set_num(m, t);
        mpz_set_ui(t, 1);
        mpz_mul_2exp(t, t, I->k + 1);
                
        mpq_set_den(m, t);
         /* fprintf(stderr, "m: "); mpq_out_str(stderr, 10, m); fprintf(stderr, " \n ");  */        
    } else {
        
        mpz_set_ui(t, 1);
        mpz_mul_2exp(t, t, -I->k - 1);
        mpz_add(t, t, I->c);
        mpq_set_num(m, t);
    }

    mpz_clear(t);
}


double slv_bintvl_get_mid_d(slv_bintvl_srcptr I)
{
    mpq_t m;
    double r;
    slv_bintvl_get_mid(m, I);
    r = mpq_get_d(m);
    mpq_clear(m);
    return r;
                     
}

void slv_add_bintvl( slv_bintvl_t* roots, 
                       slv_bintvl_srcptr I,
                       slv_info_ptr info )
{
    int b;
    int j;
    /* REMOVE THIS EVENTUALLY */
    int k;
    mpz_t c;    

    
    mpz_init_set(c, I->c);
    k = I->k;
    j = info->nb_roots;  /** The current root */
    b = info->bd; 
	
    /* debug("#roots: %lu", info->nb_roots); */
    /* debug("\t sign: %d  \t b: %d  \t k: %d", sign, b, k); */
    /* printf("***I = "); mpz_out_str(stdout, 10, c); printf("/2^%d \n", k); */

    if ( I->is_exact ) { debug("Add an exact root..."); }
    
    roots[j]->sign = info->sign;
    mpz_init(roots[j]->c);
    if (k <= b)	{
        if ( info->sign == -1 ) {
            mpz_neg(roots[j]->c, c);
            mpz_sub_ui(roots[j]->c, roots[j]->c, 1);
            mpz_mul_2exp(roots[j]->c, roots[j]->c, b-k);
        } else {
            mpz_mul_2exp(roots[j]->c, c, b-k); 
        }
        roots[j]->k = k - b;
        roots[j]->is_exact = I->is_exact;
        mpz_clear(c);
        return;
    } else {
        if ( info->sign == -1 ) {
            mpz_neg(roots[j]->c, c);
            mpz_sub_ui(roots[j]->c, roots[j]->c, 1);
        } else {
            mpz_set(roots[j]->c, c);
        }
        roots[j]->k = k - b;
        roots[j]->is_exact = I->is_exact;
    }
    mpz_clear(c);
    return;
}



int slv_is_zero_a_root(mpz_t *P, 
                       slv_bintvl_t* vec_bintvl, 
                       slv_info_ptr info)
{
    int i, j;
    slv_bintvl_t I;
    i = 0;
    /* Check if 0 is a root */
    if ( mpz_cmp_ui(P[0], 0) == 0 ) {
        i = 1; while ( mpz_cmp_ui(P[i], 0) == 0 ) { i++; }
        slv_bintvl_init(I);
        I->is_exact = 1;
        for ( j = 0; j < i; j++ ) {
            slv_add_bintvl(vec_bintvl, I, info);
            info->nb_roots++;
        }
        /* the poly has a smaller degree */
        info->dg -= i; 
        for ( j = 0; j <= info->dg; j++, i++ ) {
            mpz_set(P[j], P[i]);
        }
        slv_bintvl_clear(I);
    }
    /* debug("dg= %ld ", info->dg); */
    /* poly_pretty_print(stdout, P, info->dg); */
    return (i > 0 ? 1 : 0);
}


double get_cpu_time ()
{
    return (double)clock() /(double) CLOCKS_PER_SEC;
}

