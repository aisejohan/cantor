#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>



void set_seed(unsigned int zaadje)
{
        int fd, uit;
        unsigned int willekeurig;
        if (zaadje) {
                srand(zaadje);
                return;
        }
        fd = open("/dev/urandom", O_RDONLY);
        if (fd < 0) {
                printf("Unable to open /dev/urandom. So seed=666.\n");
                willekeurig = 666;
        } else {
                uit = read(fd, &willekeurig, sizeof(willekeurig));
                if (uit <= 0) {
                        printf("Failure reading /dev/urandom. So Seed=666.\n");
                        willekeurig = 666;
                }
        }
        srand(willekeurig);
        uit = close(fd);
        return;
}

/* Merge two sorted lists of integers */
int *merge(int *v, int *w, int flag)
{
	int a,b,i,j,k;
	int *u;

	if ((flag != 0) && (flag != 1)) {
		printf("Have to use flag=1 or 0 in merge function!\n");
		exit(1);
	}

	a = v[0];
	b = w[0];

	u = (int *) malloc((a + b + 1)*sizeof(int));
	u[0] = -1;
	i = 1;
	j = 1;
	k = 1;
	while ((i <= a) && (j<= b)) {
		if (v[i] < w[j]) {
			if (u[k-1] < v[i] + flag) {
				u[k] = v[i];
				k++;
			}
			i++;
		} else {
			if (u[k-1] < w[j] + flag) {
				u[k] = w[j];
				k++;
			}
			j++;
		}
	}
	while (i <= a) {
		if (u[k-1] < v[i] + flag) {
			u[k] = v[i];
			k++;
		}
		i++;
	}
	while (j <= b) {
		if(u[k-1] < w[j] + flag) {
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
