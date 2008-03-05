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
#include "list_degrees.h"

void test_scalars(void )
{
	scalar ma,mb,mc,md,me;

	printf("/* The prime is: %d. */\n",prime);

	ma = 0;
	printf("The number 0 is: ");
	print_scalar(ma);
	printf(".\n");
	mb = 1;
	printf("The number 1 is: ");
	print_scalar(mb);
	printf(".\n");
	mc = 541 % prime;
	printf("The number 541 is: ");
	print_scalar(mc);
	printf(".\n");
	md = sc_neg_inv(mc);
	printf("The inverse of 541 is: ");
	print_scalar(md);
	printf(".\n");
	me = (mc * md) % prime;
	printf("The number -1 is: ");
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

	random_pol(A, 100);
	random_pol(B, 100);
	random_pol(C, 100);
	pol_add(D, A, B);
	pol_mult(D, C, D);
	pol_mult(B, B, C);
	pol_mult(A, A, C);
	pol_add(C, A, B);
	times_int(C, -1, C);
	pol_add(D, C, D);

	printf("The following should be zero:\n");
	print_pol(D);
	printf("\n");

	free_pol(&A);
	free_pol(&B);
	free_pol(&C);
	free_pol(&D);
}

void test_associative_law()
{
	polynomial A,B,C,D;

	make_pol(&A);
	make_pol(&B);
	make_pol(&C);
	make_pol(&D);

	random_pol(A, 400);
	random_pol(B, 600);
	random_pol(C, 140);
	pol_mult(D, A, B);
	pol_mult(D, D, C);
	pol_mult(B, B, C);
	pol_mult(C, A, B);
	times_int(C, -1, C);
	pol_add(D, C, D);

	printf("The following should be zero:\n");
	print_pol(D);
	printf("\n");

	free_pol(&A);
	free_pol(&B);
	free_pol(&C);
	free_pol(&D);
}

void test_reduction()
{
	polynomial A,B,C,D,q,r;

	make_pol(&A);
	make_pol(&B);
	make_pol(&C);
	make_pol(&D);
	make_pol(&q);
	make_pol(&r);

	random_pol(A, 100);
	random_pol(B, 500);
	qr_reduce(r, B, q, A);
	r_reduce(C, B, A);
	times_int(C, -1, C);
	pol_add(C, r, C);
	printf("The following should be zero:\n");
	print_pol(C);
	printf("\n");
	pol_mult(C, q, A);
	pol_add(D, B, C);
	times_int(r, -1, r);
	pol_add(D, D, r);
	printf("The following should be zero:\n");
	print_pol(D);
	printf("\n");

	free_pol(&A);
	free_pol(&B);
	free_pol(&C);
	free_pol(&D);
	free_pol(&q);
	free_pol(&r);
}

void test_gcd()
{
	scalar c;
	polynomial A,B,C,D,g;

	make_pol(&A);
	make_pol(&B);
	make_pol(&C);
	make_pol(&D);
	make_pol(&g);

	random_pol(A, 20);
	random_pol(B, 20);
	random_pol(C, 20);

	gcd(D, A, B);
	pol_mult(D, D, C);

	pol_mult(A, A, C);
	pol_mult(B, B, C);
	gcd(g, A, B);
/*
	printf("The gcd is:\n");
	print_pol(g);
	printf("\n");
*/
	if (D->degree != g->degree) {
		printf("test_gcd failed\n");
		exit(1);
	}
	c = sc_neg_inv(g->coeffs[g->degree]);
	c = (D->coeffs[g->degree] * c) % prime;
	times_scalar(g, c, 0, g);
	pol_add(D, g, D);
	printf("The following should be zero:\n");
	print_pol(D);
	printf("\n");

	free_pol(&A);
	free_pol(&B);
	free_pol(&C);
	free_pol(&D);
	free_pol(&g);
}

void test_deriv()
{
	polynomial A,Ap,B,Bp,C,D;

	make_pol(&A);
	make_pol(&Ap);
	make_pol(&B);
	make_pol(&Bp);
	make_pol(&C);
	make_pol(&D);

	random_pol(A, 200);
	random_pol(B, 200);

	pol_mult(C, A, B);
	deriv(D, C);
	deriv(Ap, A);

/*
	printf("The polynomial A:\n");
	print_pol(A);
	printf("The derivative of A:\n");
	print_pol(Ap);
	printf("\n");
*/

	deriv(Bp, B);
	pol_mult(C, Ap, B);
	pol_mult(A, A, Bp);
	pol_add(C, C, A);
	times_int(C, -1, C);
	pol_add(D, C, D);

	printf("The following should be zero:\n");
	print_pol(D);
	printf("\n");

	free_pol(&A);
	free_pol(&Ap);
	free_pol(&B);
	free_pol(&Bp);
	free_pol(&C);
	free_pol(&D);
}

void test_p_power()
{
	int i,d;
	polynomial A,Ap,B,C,g;

	make_pol(&A);
	make_pol(&Ap);
	make_pol(&B);
	make_pol(&C);
	make_pol(&g);

	/* Do not make this too large. */
	d = 10;
	random_pol(A, d);

	deriv(Ap, A);
	gcd(g, Ap, A);
	if (g->degree > 0) qr_reduce(B, A, A, g);

	d = A->degree;
	if (d == 0) goto uit;
	random_pol(B, d-1);

	prime_power(C, B, A);

	printf("If it hangs at this spot then prime_power is defective!\n");
	i = 1;
	while (!equal(B, C)) {
		prime_power(C, C, A);
		i++;
	}

uit:
	free_pol(&A);
	free_pol(&Ap);
	free_pol(&B);
	free_pol(&C);
	free_pol(&g);
}

void test_print_degrees()
{
	polynomial A;

	make_pol(&A);
	random_pol(A, 1000);
	print_pol(A);
	print_degrees(A);
}


/*
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
/*
	change_prime(17);
	test_scalars();
	test_distributive_law();
	test_associative_law();
	test_reduction();
	test_gcd();
	test_deriv();
	test_p_power();
	test_print_degrees();
*/
	change_prime(547);
	test_print_degrees();
	return(0);
}
