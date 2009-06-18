/*
 *	test.c
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
#include <unistd.h>

#include "data.h"
#include "scalar.h"
#include "pol.h"
#include "xu_and_sparse.h"
#include "list_degrees.h"
#include "utils.h"

int main()
{
	int nr, i;
	int primes[100] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541};
	int *list;
	xu_polynomial A;
	sparse_polynomial sB;
	polynomial B;

	make_pol(&B);
	set_seed(0);
	random_xu(&A, 3, 3);
	printf("f = ");
	print_xu_pol(A);

	i = 0;
	while (i <= 99) {
		change_prime(primes[i]);
		nr = xu_to_sparse(&sB, A);
		if (nr) {
			sparse_to_pol(B, sB);
			printf("p = %d: ", primes[i]);
/*			printf("[");
			print_degrees_sparse(B, sB);
			printf("];\n");
*/
			list = list_degrees(B);
			print_list(list);
			free(list);
			free_sparse_pol(&sB);
		}
		i++;
	}

	free_pol(&B);
	free_xu_pol(&A);
	return(0);
}
