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
#include "xu_and_sparse.h"
#include "utils.h"
#include "emgcd.h"

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
	int i;
	int *list;
	polynomial A, B;

	make_pol(&A);
	make_pol(&B);
	random_pol(A, 5);
	copy_pol(B, A);
	i = 1;
	while (i < prime) {
		pol_mult(B, B, A);
		i++;
	}
	A->degree = 4;
	A->coeffs[0] = 0;
	A->coeffs[1] = 0;
	A->coeffs[2] = 0;
	A->coeffs[3] = 0;
	A->coeffs[4] = 1;
	pol_mult(B, A, B);
	print_pol(B);
	printf("\n");
	printf("\n");
	list = list_degrees(B);
	print_list(list);
	printf("THE END.\n");
}

void test_sparse_reduce()
{
	polynomial A, B, C, D;
	sparse_polynomial sA;

	make_pol(&A);
	make_pol(&B);
	make_pol(&C);
	make_pol(&D);
	random_pol(A, 20);
	random_pol(B, 30);

	convert_to_sparse(&sA, A);

	r_reduce(D, B, A);
	r_reduce_sparse(C, B, sA);
	times_int(C, -1, C);
	pol_add(C, C, D);
	printf("The following should be zero:\n");
	print_pol(C);

	free_pol(&A);
	free_pol(&B);
	free_pol(&C);
	free_pol(&D);
	free_sparse_pol(&sA);
}

void test_xu_to_sparse()
{
	int nr, i;
	int primes[100] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541};
	xu_polynomial A;
	sparse_polynomial sB;
	polynomial B;

	make_pol(&B);
	random_xu(&A, 4, 3);

	i = 0;
	while (i <= 10) {
		change_prime(primes[i]);
		nr = xu_to_sparse(&sB, A);
		if (nr) {
			sparse_to_pol(B, sB);
			printf("\n");
			printf("\n");
			print_pol(B);
			print_degrees_sparse(B, sB);
			free_sparse_pol(&sB);
		}
		i++;
	}

	free_pol(&B);
	free_xu_pol(&A);
}

void test_emgcd()
{
	int i;
	polynomial B, C;
	polynomial A[2];
	polynomial *EMGCD;

	make_pol(&A[0]);
	make_pol(&A[1]);
	random_pol(A[0], 10000);
	print_pol(A[0]);
	random_pol(A[1], 9999);
	print_pol(A[1]);

	EMGCD = do_emgcd(A);

	make_pol(&B);
	make_pol(&C);
	pol_mult(B, EMGCD[2], A[0]);
	pol_mult(C, EMGCD[4], A[1]);
	pol_add(B, B, C);
	times_int(B, -1, B);
	pol_add(B, B, EMGCD[0]);
	print_pol(B);
	pol_mult(B, EMGCD[3], A[0]);
	pol_mult(C, EMGCD[5], A[1]);
	pol_add(B, B, C);
	times_int(B, -1, B);
	pol_add(B, B, EMGCD[1]);
	print_pol(B);
}

int main(void )
{
	change_prime(17);
	test_scalars();
	test_distributive_law();
	test_associative_law();
	test_reduction();
	test_gcd();
	test_deriv();
	test_p_power();
	test_print_degrees();

	change_prime(5);
	test_emgcd();

/*	change_prime(61);
	test_print_degrees();
	test_sparse_reduce();
	test_xu_to_sparse(); */
	return(0);
}
