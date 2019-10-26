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

#include "info.h"
#include "utils.h"



void slv_info_init(slv_info_ptr info)
{
	info->max_depth = 0;
	info->nb_nodes = 0;
	info->nb_homo = 0;
	info->nb_trans = 0;
	info->nb_half_opt = 0;
    info->nb_pos_hack_1 = 0;
    info->nb_pos_hack_2 = 0;
    
	info->t_dg = 0;
    info->dg = 0;
	info->bd = -1;
    info->nb_roots = 0;

	return;
}


void slv_info_print(slv_info_srcptr info)
{
    printf("/*------------------------------------------------------------*/\n");
    printf("  Statistics: \n");

    printf( "  #degree : %3lu \n", info->dg);
    printf( "  #bound  : %3lu \n", info->bd);
    
    printf( "  #roots  : %3lu \n", info->nb_roots);
    printf( "  #nodes  : %3lu \n", info->nb_nodes);
    printf( "  #depth  : %3lu \n", info->max_depth);

    printf( "  #trans  : %3lu \n", info->nb_trans);
    printf( "  #homo   : %3lu \n\n", info->nb_homo);

    printf( "  #pos_h_1  : %3lu \n", info->nb_pos_hack_1);
    printf( "  #pos_h_2  : %3lu \n", info->nb_pos_hack_2);
    printf( "  #half_h   : %3lu \n", info->nb_half_opt);
    
    printf("/*------------------------------------------------------------*/\n");

}


void slv_info_print2 (slv_info_srcptr info1, slv_info_srcptr info2)
{
    printf("/*------------------------------------------------------------*/\n");
    printf("  Statistics: \n");

    printf( "   #degree : %3lu \n", info1->dg);

    printf( "             Neg\t Pos  \t Total\n");
    printf( "   #bound  : %3ld \t %3ld \n", info1->bd, info2->bd);
    printf( "   #roots  : %3lu \t %3lu  \t %3lu \n\n", info1->nb_roots, info2->nb_roots, info1->nb_roots+info2->nb_roots);

    printf( "   #nodes  : %3lu \t %3lu  \t %3lu \n", info1->nb_nodes, info2->nb_nodes, info1->nb_nodes+info2->nb_nodes);
    printf( "   #depth  : %3lu \t %3lu  \t %3lu \n", info1->max_depth, info2->max_depth, info1->max_depth+info2->max_depth);

    printf( "   #trans  : %3lu \t %3lu  \t %3lu \n", info1->nb_trans, info2->nb_trans, info1->nb_trans+info2->nb_trans);
    printf( "   #homo   : %3lu \t %3lu  \t %3lu\n\n", info1->nb_homo, info2->nb_homo, info1->nb_homo+info2->nb_homo);

    printf( "   #pos_h_1  : %3lu \t %3lu \t %3lu \n", info1->nb_pos_hack_1, info2->nb_pos_hack_1, info1->nb_pos_hack_1+info2->nb_pos_hack_1);
    printf( "   #pos_h_2  : %3lu \t %3lu  \t %3lu \n", info1->nb_pos_hack_2, info2->nb_pos_hack_2, info1->nb_pos_hack_2+info2->nb_pos_hack_2);
    printf( "   #half_h   : %3lu \t %3lu  \t %3lu \n", info1->nb_half_opt, info2->nb_half_opt, info1->nb_half_opt+info2->nb_half_opt);

    printf("/*------------------------------------------------------------*/\n");

}


void slv_info_add(slv_info_ptr info, slv_info_srcptr info1, slv_info_srcptr info2)
{
	info->dg = slv_MAX(info1->dg, info2->dg);
	info->bd = slv_MAX(info1->bd, info2->bd);
	info->nb_roots = info1->nb_roots + info2->nb_roots;
	info->nb_nodes = info1->nb_nodes + info2->nb_nodes;
	info->max_depth = slv_MAX(info1->max_depth, info2->max_depth);
	info->nb_homo = info1->nb_homo + info2->nb_homo;
	info->nb_trans = info1->nb_trans + info2->nb_trans;
	info->nb_half_opt = info1->nb_half_opt + info2->nb_half_opt;
	info->nb_pos_hack_1 = info1->nb_pos_hack_1 + info2->nb_pos_hack_1;
	info->nb_pos_hack_2 = info1->nb_pos_hack_2 + info2->nb_pos_hack_2;
}
