#include <stdlib.h>
#include <stdio.h>

/* Merge two sorted lists of integers */
int *merge(int *v, int *w)
{
	int a,b,i,j,k;
	int *u;

	a = v[0];
	b = w[0];

	u = (int *) malloc((a + b + 1)*sizeof(int));
	u[0] = -1;
	i = 1;
	j = 1;
	k = 1;
	while ((i <= a) && (j<= b)) {
		if (v[i] < w[j]) {
			if (u[k-1] <= v[i]) {
				u[k] = v[i];
				k++;
			}
			i++;
		} else {
			if (u[k-1] <= w[j]) {
				u[k] = w[j];
				k++;
			}
			j++;
		}
	}
	while (i <= a) {
		if (u[k-1] <= v[i]) {
			u[k] = v[i];
			k++;
		}
		i++;
	}
	while (j <= b) {
		if(u[k-1] <= w[j]) {
			u[k] = w[j];
			k++;
		}
		j++;
	}
	u[0] = k-1;
	return(u);
}

void print_list(int *list)
{
	int i;

	i = 1;
	printf("[");
	while (i < list[0]) {
		printf("%d, ", list[i]);
		i++;
	}
	if (list[0] > 0) {
		printf("%d]\n", list[i]);
	} else {
		printf("]\n");
	}
}
