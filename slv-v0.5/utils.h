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

#ifndef _SLV_UTILS_H_
#define _SLV_UTILS_H_


#define ALL_COEFFS_ARE_POSITIVE   -1

/* Keeps the depth of the tree. */
#define SLV_MAX_DEPTH 2000



#define get_first_1(x) 	mpz_scan1((x), 0)
#define mpz_bitsize(a) 	mpz_sizeinbase((a),2)


#define slv_MIN(a, b) ( ((a) < (b)) ? (a) : (b) )
#define slv_MAX(a, b) ( ((a) > (b)) ? (a) : (b) )

#define mpfi_sgn(I)                             \
    (mpfi_is_strictly_neg(I) ? -1 :             \
                             ( mpfi_is_strictly_pos(I) ? 1 : \
                               (mpfi_is_zero(I) ? 0 : -2)    \
                               )                             \
     )                                                       \
    
#define ucoeff(poly, n) ((poly)->coeff[n])

    
int b_verb ;

#endif /* _SLV_UTILS_H_ */
