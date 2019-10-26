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


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#include "vca_solver_1.h"
#include "utils.h"


/*
  Descartes' test based on Taylor shift and number of sign variations.
  Basic routine with no optimizations.
  Just for debugging.
 */
unsigned long Descartes_test_3( mpz_t* P, unsigned long deg, long sigh, long *flag)
{
	mpz_t* Q;
	long i;
    long v; 

	Q = (mpz_t*) malloc( (deg + 1) * sizeof( mpz_t));
	for (i = 0; i <= deg; i++) {
	    mpz_init_set( Q[deg-i], P[i]);
	}

	slv_poly_taylor_shift_by_1( Q, deg);
    v = var( Q, deg);

	for (i = 0; i <= deg; i++)
		mpz_clear(Q[i]);
	free(Q);

	return v;
}


/* 
   Upper bound on the number of positive real roots using Descartes' rule of sign/
   Optimized procedure using various hacks described in [Rouillier,Zimmermann:JCAM:2004]
   This is a slightly modified version of a varariant of 
   [Hanrot, Rouillier, Zimmermann, and Petitjean].
*/
inline
int
Descartes_test(mpz_t *P, unsigned long deg, long sigh,
               long* status,
               slv_info_ptr info,
               mpz_t* Q)
{
    
	unsigned long V = 0;
	long i, j, s, t;
    
    info->nb_nodes++;
	j = deg; t = mpz_sgn(P[deg]);

    /* poly_print_list(stdout, P, deg); */
    while ( (j >= 0) && ( (mpz_sgn(P[j]) == t) || 
                          (mpz_sgn(P[j]) == 0) ) )  {
        j--;
	}
	if ( j < 0 ) {
        info->nb_pos_hack_1++;
        *status = ALL_COEFFS_ARE_POSITIVE;
        return V;
	}

    SLV_POLY_SET(Q, P, deg, i);	
	for ( j = 0; j <= deg-1; j++ ) { mpz_add( Q[j+1], Q[j+1], Q[j]); }

    s = mpz_sgn(Q[deg]);
    
    *status = s && (s == mpz_sgn( P[0])) && (s == -sigh);
    
	for ( i = 1; i <= deg-1; i++ ) {
        j = deg - i; t = s;
        while ( (j >= 0) && (t == 0) ) { t = mpz_sgn(Q[j]); j--; }
        while ( (j >= 0) && ( (mpz_sgn(Q[j]) == t) ||
                              (mpz_sgn(Q[j]) == 0) ) ) {  j--;  }
		if ( j < 0 ) {
            info->nb_pos_hack_2++;
            return V;
		}
        
		for ( j = 0; j <= deg - i - 1; j++ ) { mpz_add(Q[j+1], Q[j+1], Q[j]); }
        
		if ( s == 0 ) {
			s = mpz_sgn(Q[deg-i]);
		} else {
			if ( s == -mpz_sgn(Q[deg-i]) ) {
				if ( ((V == 1) && !*status) || (V == 2) ) {
                    return (V + 1);
				}
				V++; s = -s;
			}
		}
	}
	if ( s == -mpz_sgn(Q[0]) ) V++;
    
	return V;
}

 


/* Memory efficient variant. Perform all the computations with 2 polynomials.  */

void VCA_0_1_mem(mpz_t *FF, 
				unsigned long dg,
                 slv_bintvl_t* vec_bintvl, 
                 slv_info_ptr info)
{

    long i, V;
    long k = 0;
    
    mpz_t c;
    
    long status = 0;
    long shalf= 1;
    
    info->dg = dg;

    long cnt = 0;
    
    mpz_t* P;
    mpz_t* Q;
    slv_lst_bintvl  queue;
   
    slv_bintvl_ptr I;
    I = malloc(sizeof(slv_bintvl_t)); 

    SLV_POLY_INIT_SET(P, FF, dg, i);
    SLV_POLY_INIT(Q, dg, i);
    
    // V = Descartes_test(P, dg, shalf, &status, info, Q);
    V = var(P, dg);
    /* debug( "var: %ld \n", V); */
      
    if ( V == 0 ) { return ; }

    slv_bintvl_init(I);
    if ( V == 1 ) {
        slv_add_bintvl(vec_bintvl, I, info);
        info->nb_roots++;            
        return;
    }
    /* The poly has more than one sign variations. Split and continue.  */
    
    SLIST_INIT(&queue);

    
    slv_bintvl_ptr Il = malloc(sizeof(slv_bintvl_t));
    slv_bintvl_ptr Ir = malloc(sizeof(slv_bintvl_t));
    
    Il = malloc(sizeof(slv_bintvl_t));         
    slv_bintvl_init(Il);
    Il->k = 1;
    info->max_depth = 1;
    /* Il->depth = 1; */

    Ir = malloc(sizeof(slv_bintvl_t));                
    slv_bintvl_init(Ir);
    Ir->k = 1;
    mpz_set_ui(Ir->c, 1);
    /* Ir->depth = 1; */
                 
    SLIST_PUSH(queue, Ir);
    SLIST_PUSH(queue, Il);
    
    k = 0;
    mpz_init_set_ui(c, 0);

    while ( (!SLIST_EMPTY(&queue)) && info->max_depth <= SLV_MAX_DEPTH)  {
        SLIST_POP(queue, I);
        
        info->max_depth = slv_MAX(info->max_depth, I->k);
        
        /* Construct the polynomial depending on the previous k */
        if ( k < I->k ) {
            slv_poly_scale_2exp(P, dg, -1);
            info->nb_homo++;
            k++;
        } else if ( k == I->k) {
            slv_poly_taylor_shift_by_1(P, dg);
            info->nb_trans++;
        } else if ( k > I->k ) {
            slv_poly_taylor_shift_by_1(P, dg);
            info->nb_trans++;
            slv_poly_scale_2exp(P, dg, k - (I->k));  
            info->nb_homo++;
            k = I->k;
        }

        /*
          Compute the sign of P(1/2) ; thus if Descartes bound is 2,
          whereas sign(P(0)) = sign(P(1)) = -sign(P(1/2)) we have
          found two roots.
        */
        shalf = slv_poly_sgn_eval_at_half(P, dg); 
            
        V = Descartes_test(P, dg, shalf, &status, info, Q);
        switch ( V ) {
        case 0:
            slv_bintvl_clear(I);
            free(I);
          break;
                
        case 1:
          slv_add_bintvl(vec_bintvl, I, info);
          info->nb_roots++;            
          
          slv_bintvl_clear(I);
          free(I);
          break;
                
        case 2:
            if (status) {
                /* debug("status: %ld", status); */
                       
                /* There is a root in (0, 1/2) and (1/2, 1) */
                mpz_mul_2exp(I->c, I->c, 1);
                I->k++;
                slv_add_bintvl(vec_bintvl, I, info);
                info->nb_roots++;
                
                mpz_add_ui(I->c, I->c, 1);
                slv_add_bintvl(vec_bintvl, I, info);
                info->nb_roots++;

                info->nb_half_opt++;

                slv_bintvl_clear(I);
                free(I);
                break;
            }
                
        default:
            /* The left child */
            Il = malloc(sizeof(slv_bintvl_t));
            slv_bintvl_init_set(Il, I);
            Il->k = I->k + 1;
            mpz_mul_2exp(Il->c, Il->c, 1);
                 
            /* The right child */
            Ir = malloc(sizeof(slv_bintvl_t));                
            slv_bintvl_init_set(Ir, I);
            Ir->k = I->k + 1;
            mpz_set(Ir->c, Il->c);
            mpz_add_ui(Ir->c, Ir->c, 1);

            SLIST_PUSH(queue, Ir);
            SLIST_PUSH(queue, Il);
        }
        cnt++; /* For debug purposes only */
    }
    check_debug(SLIST_EMPTY(&queue), "There are intervals that we DID NOT consider!");

    mpz_clear(c);
    SLV_POLY_CLEAR(Q, dg, i);
    SLV_POLY_CLEAR(P, dg, i);
    return;

 SLV_ERROR exit(-1);
}

 

 
/* 
   TODO: We should really use a structure for the statistics.
*/

slv_bintvl_t*
VCA(mpz_t* A, unsigned long dg, slv_info_ptr info, int algo)
{

	slv_bintvl_t* vec_bintvl = (slv_bintvl_t*) malloc(dg * sizeof(slv_bintvl_t));

	int j;
	mpz_t *F;

    info->dg = dg;
    /* TODO:
       - Isolate the POS and NEG roots.
       - Check for symmetric coefficients (only the pos roots)
       - Option for the power hack (in this case we need a fast refinement routine).
     */

    SLV_POLY_INIT_SET(F, A, dg, j);

    /* Is 0 a root ? */
    slv_is_zero_a_root(F, vec_bintvl, info);
    /* debug("dg= %ld ", dg); */
    /* poly_pretty_print(stdout, F,dg); */
    

	long k = slv_poly_root_upper_bound_2exp(F, dg);
	slv_poly_scale_2exp( F, dg, k);

    info->bd = k;

	/* printf( "b = 2^%ld \n", info->bd); */
    /* poly_pretty_print(stdout, F, dg); */

    /* printf("//---- START NEW SOLVER -------------------//\n\n\n"); */
    if (algo == 1) {
        VCA_0_1_mem(F, dg, vec_bintvl, info);
    } else {
        fprintf(stderr, "No such algo !!! \n");
        exit(-1);
    }
    /* debug("#roots: %lu", info->nb_roots); */
    /* printf("//---- END NEW SOLVER  --------------------//\n\n\n\n"); */
   
    /* poly_pretty_print(stdout, F, dg); */
    
    SLV_POLY_CLEAR(F, dg, j); 

	printf( "Finish isolating \n");
	return vec_bintvl;
}






void print_root(FILE *stream, slv_bintvl_t z)
{
    mpz_t tmp;
    mpz_init( tmp);


    mpz_t c1, c2;
    long k1 = 0;
    long k2 = 0;

    mpz_init(c1);
    mpz_init(c2);

    mpz_init_set(c1, z->c);
    
    if (z->k <= 0){ 
        /* mpz_out_str(stream, 10, z->c);  */
    } else {
        k1 = z->k;
    }
    
    if (z->k <= 0) {
        mpz_set_ui(c2, 1);
        mpz_mul_2exp(c2, c2, -z->k);
        mpz_add(c2, z->c, c2);
    } else {
        mpz_add_ui(c2, z->c, 1);
        k2 = z->k;
    }
    
    
    fprintf(stream, "(");
    mpz_out_str(stream, 10, c1); fprintf(stream, "/2^%ld", k1);
    fprintf(stream, ", ");  
    mpz_out_str(stream, 10, c2); fprintf(stream, "/2^%ld", k2);
    fprintf(stream, ")");

    return ;

    if ( z->is_exact ) {}
    
    if (z->is_exact != 1) { fprintf(stream, "("); }

    if (z->k <= 0){ 
        mpz_out_str(stream, 10, z->c); 
    } else { 
        mpz_out_str(stream, 10, z->c); 
        fprintf(stream, "/2^%ld", z->k); 
    }

    if (z->is_exact == 1) { printf("\n"); return; }

    fprintf(stream, ", ");  


    if (z->k <= 0) {
        mpz_set_ui(tmp, 1);
        mpz_mul_2exp(tmp, tmp, -z->k);
        mpz_add(tmp, z->c, tmp);
        mpz_out_str(stream, 10, tmp);
    } else {
        mpz_add_ui(tmp, z->c, 1);
        mpz_out_str(stream, 10, tmp); fprintf(stream, "/2^%ld", z->k);
    }

    fprintf(stream, ") \n");
  
    mpz_clear(tmp);
}



void print_roots_all(FILE* stream, slv_bintvl_t* roots, long nbr)
{
	fprintf(stream, "#roots = %lu \n", nbr);
	fprintf(stream, "[\n");
	for (long i = 0; i < nbr; i++) {
		fprintf(stream, " %2d   ", roots[i]->sgn_left);
		print_root(stream, roots[i]);
		fprintf(stream, "\n");
	}
	fprintf(stream, "]\n");
}




void slv_adjust_signs(mpz_t* P, long dg, slv_bintvl_t* roots, slv_info_srcptr info)
{
    int s;
    long i;

    
    for (i = 0; i < info->nb_roots; i++) {
        if ( roots[i]->is_exact ) { continue; }
        
        if ( roots[i]->k > 0 ) {
            s = slv_poly_sgn_eval_at_c_2exp(P, dg, roots[i]->c, roots[i]->k);
        } else {
            s = slv_poly_sgn_eval_at_c(P, dg, roots[i]->c);
        }
        
        roots[i]->sgn_left = s;
    }
}



int slv_refine(mpz_t* P, long dg, slv_bintvl_t I, long t)
{
    int sgn_m;
    long i;
    
    mpz_t m;
    mpz_init(m);

    /* fprintf(stderr, "c: ");mpz_out_str(stderr, 10, I->c); fprintf(stderr, " ");  */
    // fprintf(stderr, "2^%ld \n", I->k);
    for (i=0; i < t; i++) {
       /* fprintf(stderr, "bm: "); mpz_out_str(stderr, 10, I->c); fprintf(stderr, " "); */
       /* fprintf(stderr, " b  2^%ld \n", I->k); */
        if (I->k < 0) {
            mpz_set_ui(m, 1);
            mpz_mul_2exp(m, m, -I->k - 1);
            mpz_add(m, m, I->c);
        } else {
            mpz_mul_ui(m, I->c, 2);
            mpz_add_ui(m, m, 1);
        }
        
        /* fprintf(stderr, "m: "); mpz_out_str(stderr, 10, m); fprintf(stderr, " "); */
        /* fprintf(stderr, "2^%ld \n", I->k); */
        sgn_m = slv_poly_sgn_eval_at_c_2exp(P, dg, m, I->k+1);
        
        
        /* debug("\n\t Signs:  %d %d  %d",  I->sgn_left, sgn_m,  -I->sgn_left); */
        if (sgn_m == 0) {
            mpz_set(I->c, m);
            I->k++;
            I->is_exact = 1;
            mpz_clear(m);
            return 0;
        }

        if (sgn_m == I->sgn_left) {
            mpz_set(I->c, m);
            I->k++;
        } else if (sgn_m == -I->sgn_left) {
            if (I->k >= 0) { mpz_mul_ui(I->c, I->c, 2);}
            I->k++;
        } else {
            debug("We should never be here!  %d %d  %d",  I->sgn_left, sgn_m,  -I->sgn_left);
            return -1;
        }
           
    }
    mpz_clear(m);
    return 1;
}


int slv_refine_until(mpz_t* P, long dg, slv_bintvl_t I, long t)
{
    // t should be negative
    // debug("here");
    check_debug(t >= 0, "The width should be small");
    long m;

    if (t >= I->k) {
        m = t - I->k;
        /* debug("m : %ld \t %ld", m, I->k); */
        return slv_refine(P, dg, I, m);
    }
    return 1;
    
    SLV_ERROR exit(-1);
}



void slv_refine_until_all(mpz_t* P, long dg, slv_bintvl_t* roots, long nbr, long t)
{
  long i;
	for (i = 0; i < nbr; i++) {
		slv_refine_until(P, dg, roots[i], t);
	}
}




slv_bintvl_t* VCA_0_infinity(mpz_t* A, unsigned long dg, slv_info_ptr info, int algo)
{
	int j;
	mpz_t *F;

	info->dg = dg;
	slv_bintvl_t* roots = (slv_bintvl_t*) malloc(dg * sizeof(slv_bintvl_t));

    SLV_POLY_INIT_SET(F, A, dg, j);

    check_debug(mpz_sgn(A[0]) != 0, "You should have already checked for 0 as a root!");
    /* slv_is_zero_a_root(F, roots, info_pos); */

	long k = slv_poly_root_upper_bound_2exp(F, dg);
    info->bd = k;
	slv_poly_scale_2exp(F, dg, k);
  

    if (algo == 1) {
        VCA_0_1_mem(F, dg, roots, info);
    } else {
        fprintf(stderr, "No such algo !!! \n");
        exit(-1);
    }
    
    SLV_POLY_CLEAR(F, dg, j);
    slv_adjust_signs(A, dg, roots, info);

	//printf( "Finish isolating \n");
	return roots;

    SLV_ERROR exit(-1);
}

 

slv_bintvl_t*
VCA_all(mpz_t* A, unsigned long dg, slv_info_ptr info_pos, slv_info_ptr info_neg, int algo)
{
	int j;
	mpz_t *F;

	info_pos->dg = dg;


	slv_bintvl_t* roots_pos = (slv_bintvl_t*) malloc(dg * sizeof(slv_bintvl_t));

    SLV_POLY_INIT_SET(F, A, dg, j);

    /* Is 0 a root ? */
    slv_is_zero_a_root(F, roots_pos, info_pos);


    // Isolate the positive roots 
    info_pos->sign = 1;
    roots_pos = VCA_0_infinity(F, dg, info_pos, algo);
    slv_adjust_signs(A, dg, roots_pos, info_pos);

    // debug("Finish with POS. # = %ld", info_pos->nb_roots);



    // Isolate the negative roots
	slv_bintvl_t* roots_neg = (slv_bintvl_t*) malloc(dg * sizeof(slv_bintvl_t));
    
    SLV_POLY_INIT_SET(F, A, dg, j);
    for (j = 1; j <= dg; j += 2) { mpz_neg(F[j], F[j]); }
    info_neg->sign = -1;
    roots_neg = VCA_0_infinity(F, dg, info_neg, algo);
    slv_adjust_signs(A, dg, roots_neg, info_neg);      
    
    // debug("Finish with NEG. # = %ld", info_neg->nb_roots);
   
    // debug( "roots %ld + %ld ", info_pos->nb_roots, info_neg->nb_roots);

 
    SLV_POLY_CLEAR(F, dg, j);


    // Merge the results 

    slv_bintvl_t* roots = (slv_bintvl_t*) malloc(dg * sizeof(slv_bintvl_t));

    for (j=0; j< info_neg->nb_roots; j++) {
        slv_bintvl_init_set(roots[j], roots_neg[info_neg->nb_roots-j-1]);
        // print_root(stdout, roots[j]);
        //  printf("\n");
    }
    for (j=0; j< info_pos->nb_roots; j++) {
        slv_bintvl_init_set(roots[info_neg->nb_roots+j], roots_pos[j]);
    }   

        
	printf( "Finish isolating \n");
	return roots;
}

