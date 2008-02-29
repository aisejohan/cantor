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

void test_scalars(void )
{
	scalar ma,mb,mc,md,me;

	change_prime(7);

	printf("/* The prime is: %d. */\n",prime);

	ma = 0;
	printf("The number 0 is: ");
	print_scalar(ma);
	printf(".\n");
	mb = 1;
	printf("The number 1 is: ");
	print_scalar(mb);
	printf(".\n");
	mc = 541;
	printf("The number 541 is: ");
	print_scalar(mc);
	printf(".\n");
	md = sc_inv(mc);
	printf("The inverse of 541 is: ");
	print_scalar(md);
	printf(".\n");
	me = mc * md % prime;
	printf("The number 1 is: ");
	print_scalar(me);
	printf(".\n");
	printf("\n");

}

void test_distributive_law(void )
{
	polynomial A,B,C,D;

	make_pol(&A);
	make_pol(&B);
	make_pol(&C);
	make_pol(&D);
printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n");

	random_pol(A,10);
	print_pol(A);
	random_pol(B,10);
	print_pol(B);
	random_pol(C,10);
	pol_add(D,A,B);
	print_pol(D);
	pol_mult(D,C,D);
	pol_mult(B,B,C);
	pol_mult(A,A,C);
	pol_add(C,A,B);
	times_int(C,-1,C);
	pol_add(D,C,D);

	printf("The following should be zero:\n");
	print_pol(D);
	printf("\n");

	free_pol(&A);
	free_pol(&B);
	free_pol(&C);
	free_pol(&D);
}


/*
void test_deriv()
{
	polynomial A,B,C,D;
	make_pol(A);
	make_pol(B);
	make_pol(C);
	make_pol(D);

	random_pol(A,10);
	random_pol(B,200);
	pol_mult(A,B,C);
	deriv(C,C);
	times_int(-1,C,C);
	deriv(A,D);
	pol_mult(B,D,D);
	pol_add(D,C,D);
	deriv(B,B);
	pol_mult(A,B,B);
	pol_add(D,B,D);

	printf("The following should be zero:\n");
	print_pol(D);
	printf("\n");

	free_pol(A);
	free_pol(B);
	free_pol(C);
	free_pol(D);
}

void test_exponents()
{
	int i;
	struct exponents e;
	polynomial A,B,C,D;

	make_pol(A);
	make_pol(B);
	make_pol(C);
	make_pol(D);
	
	random_pol(B,8);
	for(i = 1; i <= 3; i++) {
		print_pol(B);
		e = take_exponents(B);
		printf("The length is %d and the degree is %d.\n",
			B->length, B->degree);
		printf("The valuation is %d and the pivot is %d.\n",
			e.val, e.piv);
		times_int(p*p+1,B,B);
	}
	printf("\n");

	free_pol(A);
	free_pol(B);
	free_pol(C);
	free_pol(D);
}


void test_associative_law()
{
	polynomial A,B,C,D;

	make_pol(A);
	make_pol(B);
	make_pol(C);
	make_pol(D);

	random_pol(A,400);
	random_pol(B,600);
	random_pol(C,140);
	pol_mult(A,B,D);
	pol_mult(C,D,D);
	pol_mult(A,B,B);
	pol_mult(B,C,C);
	times_int(-1,C,C);
	pol_add(D,C,D);

	printf("The following should be zero:\n");
	print_pol(D);
	printf("\n");

	free_pol(A);
	free_pol(B);
	free_pol(C);
	free_pol(D);
}

void test_reduction()
{
	polynomial A,B,C,D;

	make_pol(A);
	make_pol(B);
	make_pol(C);
	make_pol(D);

	random_pol(A,400);
	random_pol(B,500);
	copy_pol(B,D);
	reduce_by_g_q(B,A,C);
	pol_mult(C,A,C);
	pol_add(C,B,C);
	times_int(-1,C,C);
	pol_add(C,D,D);
	printf("The following should be zero:\n");
	print_pol(D);
	printf("\n");

	free_pol(A);
	free_pol(B);
	free_pol(C);
	free_pol(D);
}

void test_inverses()
{
	polynomial A,B,C,D;

	make_pol(A);
	make_pol(B);
	make_pol(C);
	make_pol(D);

	random_pol(B,17);
	times_int(p,B,B);
	random_pol(A,0);
	pol_add(A,B,A);
	invert_unit(A,C);
	printf("The degree of the inverse is %d.\n", C->degree);
	pol_mult(A,C,D);
	
	printf("The following should be one:\n");
	print_pol(D);
	printf("\n");

	free_pol(A);
	free_pol(B);
	free_pol(C);
	free_pol(D);
}

void test_substitution()
{
	int i,uit;
	polynomial A,B,C,D;

	make_pol(A);
	make_pol(B);
	make_pol(C);
	make_pol(D);

	random_pol(A, 20);
	random_pol(B, 20);
	substitute(A, B, D);
	deriv(D, D);
	deriv(A, A);
	substitute(A, B, C);
	deriv(B, B);
	pol_mult(B, C, C);
	negate_pol(C, C);
	pol_add(C, D, D);

	printf("The following should be zero:\n");
	print_pol(D);
	printf("\n");

	for(i = 0; i<= 3; i++) {
		random_pol(A, 50);
		random_pol(B, 50);
		uit = euclidean_mod_p(A, B, C, D);
		if (!uit) {
			printf("Not rel prime: \n");
			print_pol(A);
			print_pol(B);
		} else {
			printf("Rel prime: \n");
			invert_unit_mod_g(A, B, C);
			print_pol(A);
			print_pol(B);
			print_pol(C);
		}
	}
	printf("\n");

	free_pol(A);
	free_pol(B);
	free_pol(C);
	free_pol(D);
}

void test_frob_lift()
{
	polynomial A,B,C,D;

	make_pol(A);
	make_pol(B);
	make_pol(C);
	make_pol(D);

	random_pol(A,20);
	print_pol(A);

	frob_lift(A, B);
	print_pol(B);
	put_out(B);
	printf("\n");

	free_pol(A);
	free_pol(B);
	free_pol(C);
	free_pol(D);
}

void test_substitute_speed()
{
	polynomial A,B,C,D;

	make_pol(A);
	make_pol(B);
	make_pol(C);
	make_pol(D);

	random_pol(A, 251);
	random_pol(B, 36);
	substitute(B, A, C);
	printf("A = ");
	print_pol(A);
	printf("B = ");
	print_pol(B);
	printf("C = ");
	print_pol(C);

	free_pol(A);
	free_pol(B);
	free_pol(C);
	free_pol(D);
}
*/

int main(void )
{
	test_scalars();
	test_distributive_law();
	return(0);
}
