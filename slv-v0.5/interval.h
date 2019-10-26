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

#ifndef INTERVAL_H_
#define INTERVAL_H_


#include <stdio.h>
#include <gmp.h>
#include <sys/types.h>
#include <sys/resource.h>


#include "dbg.h"
#include "sys_queue.h"
#include "info.h"



/** \struct slv_bintvl
    \brief Bisection Interval struct.
    
    Struct for keeping an interval of a subdivision algorithm.
 */
struct slv_bintvl
{
    mpz_t         c;        /**< The interval is \f$(c/2^k, (c+1)/2^k )\f$ */
    long          k;         

    unsigned int  is_exact;
    int           sgn_left;	    /**< sign on the left endpoint */

    unsigned int   sign; 	    /**< are we solving in a negative interval? */

    mpz_t*         P;     		/**< polynomial at the current node. */

    SLIST_ENTRY(slv_bintvl) pnext; /**< pointer to the next entry in the list. */
};



/* List of Intervals */
SLIST_HEAD(slv_lst_intvl_t, slv_bintvl) ;
typedef struct slv_lst_intvl_t slv_lst_bintvl;


/* Abbreviation for basic list operations */
#define SLIST_PUSH(queue, I)  SLIST_INSERT_HEAD(&queue, I, pnext);

#define SLIST_POP(queue, I)                       \
    I = SLIST_FIRST(&queue);                      \
    SLIST_REMOVE_HEAD(&queue, pnext);             \

         

typedef struct slv_bintvl 			slv_bintvl_t[1];
typedef struct slv_bintvl* 			slv_bintvl_ptr;
typedef const struct slv_bintvl* 	slv_bintvl_srcptr;


void slv_bintvl_print    (slv_bintvl_srcptr I);
void slv_bintvl_init     (slv_bintvl_ptr a);
void slv_bintvl_init_set (slv_bintvl_ptr a, slv_bintvl_srcptr b);
void slv_bintvl_clear    (slv_bintvl_ptr a);

void slv_bintvl_get_mid(mpq_t m, slv_bintvl_srcptr I);
double slv_bintvl_get_mid_d(slv_bintvl_srcptr I);



void slv_add_bintvl(slv_bintvl_t* vec_bintvl, 
                      slv_bintvl_srcptr I,
                      slv_info_ptr info );

int slv_is_zero_a_root(mpz_t *P, 
                       slv_bintvl_t* vec_bintvl, 
                       slv_info_ptr info);


double get_cpu_time();



#endif /* INTERVAL_H_ */
