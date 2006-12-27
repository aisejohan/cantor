/*
 *	scalar.c
 *
 * 	Copyright 2006 Johan de Jong
 *
 *	This file is part of Frobenius
 *
 *	Frobenius is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation; either version 2 of the License, or
 *	(at your option) any later version.
 *
 *	Frobenius is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with Frobenius; if not, write to the Free Software Foundation, 
 *	Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *									*/

#include <stdio.h>
#include <stdlib.h>
#include "data.h"
#include "scalar.h"

/* Only called once. */
void setup_scalars(void)
{
	mpz_init_set_ui(prime, (unsigned long) p);
	mpz_init(modulus);
	mpz_ui_pow_ui(modulus, (unsigned long) p, (unsigned long) r);
}

void printmscalar(mscalar a)
{
	mpz_cdiv_q_ui(temp,modulus,2);
	if (mpz_cmp(a, temp)>0) {
		mpz_sub(temp, a, modulus);
		mpz_out_str(stdout, (int) 10, temp);
	} else {
		mpz_out_str(stdout, (int) 10, a);
	}
}

/* Divides a by b. If b is not a unit then this assumes 	*
 * valuation(a) >= valuation(b), and the result is lifted	*
 * to an integer mod p^r.					*/
/* Does not destroy a and b.					*/
void sc_div(mscalar a, mscalar b, mscalar c)
{
	unsigned long e;
	
	e = mpz_remove(c, b, prime);
	mpz_invert(c, c, modulus);
	mpz_mul(c, c, a);
	mpz_mod(c, c, modulus);
	while (e) {
		mpz_divexact_ui(c, c, (unsigned long) p);
		e--;
	}
}
