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

#ifndef __SLV_INFO_H__
#define __SLV_INFO_H__


#include <stdio.h>
#include <gmp.h>
#include <sys/types.h>
#include <sys/resource.h>


#include "dbg.h"
#include "sys_queue.h"



struct slv_info
{
	unsigned long 	max_depth;
	unsigned long 	nb_nodes;
	unsigned long	nb_homo;
	unsigned long 	nb_trans;
	unsigned long 	nb_half_opt;
    unsigned long   nb_pos_hack_1;
    unsigned long   nb_pos_hack_2;

	unsigned long 	t_dg;       /** The degree of the input */
    unsigned long 	dg;         /** The current degree. It might not be equal 
                                    to the t_dg if we find a rational root 
                                    in the process of isolation  */
    
    int             sign;       /** -1 if we solve for negative roots  */       
	long			bd;         /** The root bound is 2^{bd} */
    unsigned long	nb_roots;   /** THe number of roots */

} ;


typedef struct slv_info 			slv_info_t[1];
typedef struct slv_info* 			slv_info_ptr;
typedef const struct slv_info* 	    slv_info_srcptr;

void slv_info_init(slv_info_ptr info);
void slv_info_print(slv_info_srcptr info);

void slv_info_print2 (slv_info_srcptr info1, slv_info_srcptr info2);

void slv_info_add(slv_info_ptr info, slv_info_srcptr info1, slv_info_srcptr info2);

#endif /* __SLV_INFO_H__ */
