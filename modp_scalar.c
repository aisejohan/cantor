/*
 *	modp_scalar.c
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
#include "modp_scalar.h"

/* Only called once. */
void setup_modp_scalars(void)
{
	unsigned int i,j;
	
	NEG[0] = 0;
	for(i = 1; i <= p-1; i++) NEG[i] = p-i;
	
	INV[0] = -1; /* Garanteed to cause trouble! */
	for(i = 1; i <= p-1; i++) {
		j=1;
		while (j*i % p != 1) j++;
		INV[i] = j;
	}

	for(i = 0; i<= p-1; i++) {
		for(j = 0; j <= p-1; j++) {
			ADD[i][j] = (i + j) % p;
			MUL[i][j] = (i*j) % p;
		}
	}
}

void print_modp_scalar(unsigned int a)
{
	printf("%u", a);
}

void test_print(void)
{
	int i,j;

	printf("And here are the results of the Swedish jury.\n");
	printf("Negatives:\n");
	for(i = 0; i <= p-1; i++) printf("%u, ", NEG[i]);
	printf("\n");
	printf("Inverses:\n");
	for(i = 0; i <= p-1; i++) printf("%u, ", INV[i]);
	printf("\n");
	printf("Additives:\n");
	for(i = 0; i <= p-1; i++) {
		for(j = 0; j <= p-1; j++) printf("%u, ", ADD[i][j]);
		printf("\n");
	}
	printf("\n");
	printf("Multiplicatives:\n");
	for(i = 0; i <= p-1; i++) {
		for(j = 0; j <= p-1; j++) printf("%u, ", MUL[i][j]);
		printf("\n");
	}
	printf("\n");
}
