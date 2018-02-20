#include<stdlib.h>
#include<stdio.h>
#include"quadform.h"

#define N 1000

double drand()
{
    int r;
    r = rand();
    return (double) r / RAND_MAX;
}    

int main(int argc, char *argv[])
{
    double A[N][N], v[N], ans;
    int i, j;

    for (i=0; i<N; i++) {
        v[i] = drand();
        for (j=0; j<N; j++)
            A[i][j] = drand();
    }
    ans = quadform(N, A, v);
    printf("Value of quadratic form = %20.14e\n", ans);
    return 0;
}
