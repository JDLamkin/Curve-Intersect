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
#include <stdbool.h>

#include "dbg.h"

#include "info.h"
#include "interval.h"

#include "poly_ops.h"
#include "vca_solver_1.h"



void usage(void)
{
	char progname[] = "test_mpzx";

	printf("\n");
    printf("SLV version 0.5 \n\n");
	printf( "The syntax is: \n");
	printf( "   %s  [-h] [-f file] [-i prec] [-p] \n", progname);
	printf( "Details:\n");
    printf( " -h       : print help message\n");
	printf( " -f fname : read the coeffs from the fname (default file is ./test/dat) \n");
	printf( " -i prec  : the output intervals have width 2^(-prec) \n");
	printf( " -p       : print roots (Default: no print) \n");
	// printf( " -v       : Verbocity level (default = 0) \n");


	printf("\n\n");
}

int main(int argc, char* argv[])
{
    mpz_t* F;
    int i;
    unsigned long d;
    FILE *file_in = stdin;    
    char dfname[] = "./test.dat";
    char* fname;
    
    double st;


    
	long ref_prec = 10;
	bool b_print_roots = false;
    fname = dfname;

	while ((argc > 1) ) {
		if (strcmp(argv[1], "-h") == 0)	{
			usage();
			return 0;
		}
		if (strcmp(argv[1], "-f") == 0)	{
			// printf("1. %s\n", argv[2]);
			fname = argv[2];
		}
		if (strcmp(argv[1], "-i") == 0)	{
			// printf("2. %s\n", argv[2]);
			ref_prec = atoi(argv[2]);
		}
		if (strcmp(argv[1], "-p") == 0)	{
			b_print_roots = true;
		}
		/* if (strcmp(argv[1], "-v") == 0)	{ */
		/* 	// printf("2. %s\n", argv[2]); */
		/* 	b_verb = atoi(argv[2]); */
		/* } */

		++argv;
		--argc;
	}

    // Check if fname is NULL


	if ((file_in = fopen(fname, "r")) == NULL)  {
		fprintf(stderr, "Unable to open file %s.\n", argv[0]);
		exit(-1);
    }
    F = slv_poly_read_from_file( file_in, 0, &d);
   /* slv_poly_pretty_print( stdout, F, d); */

    printf("SLV version 0.5 \n\n");
    printf("Isolating %s \n ", fname);
    
    /* Vector of root intervals */
    slv_bintvl_t* roots;

    /* info */
	slv_info_t info_pos, info_neg;
	slv_info_init(info_pos);
	slv_info_init(info_neg);


    /* printf( "degree = %ld\n", d); */


    st = get_cpu_time();
    roots = VCA_all(F, d, info_pos, info_neg, 1);
    st = get_cpu_time() - st;

    // print the statistics
    slv_info_print2(info_neg, info_pos);



    // debug("Refinement process...");
    /* poly_pretty_print(stderr, F, info->dg); */
    printf("Refinement process... \n");
    

	double  rf_time = 0.0;
	rf_time = get_cpu_time();
	slv_refine_until_all(F, d, roots, info_pos->nb_roots + info_neg->nb_roots, ref_prec);
	rf_time = get_cpu_time() - rf_time;

	if (b_print_roots) {
		print_roots_all(stdout, roots,  info_pos->nb_roots + info_neg->nb_roots);
	}

    fclose( file_in);
    SLV_POLY_CLEAR(F, d, i);
    
    
    printf("\n");
    printf("Solving    time: %f s\n", st);
    printf("Refinement time: %f s\n", rf_time);
    printf("Total      time: %f s\n\n", st + rf_time);
    printf("\n");
    
    return 0;
}



