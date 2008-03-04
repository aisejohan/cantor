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

static scalar *invs=NULL;

void change_prime(int p)
{
	int i,j;

	if (p <= 0) {
		printf("Primes are positive. Stop.\n");
		exit(1);
	}
	if (invs) free(invs);

	invs = (scalar *)malloc(p*sizeof(scalar));
	prime = p;
	i = 1;
	while (i < p) {
		j = 1;
		while (1 != ((i*j) % prime)) j++;
		invs[i] = j;
		i++;
	}
}

void print_scalar(scalar a)
{
	printf("%d", a);
}

scalar sc_inv(scalar i)
{
	i = i % prime;
	return(invs[i]);
}
