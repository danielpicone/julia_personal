#include<stdio.h>
#include<cblas.h>

void print_vector(int lenx, double* x, char* label)
{
    int j;
    printf("%s", label);
    for(j=0; j<lenx; j++){
        printf("%8.4f", x[j]);
    }
    printf("\n");
}

int main(int argc, char* argv[])
{
    double x[] = { 1.0, 0.0, -3.0, 2.0 };
    double y[] = { -2.0, 1.0, 5.0, 0.0 };
    double a = 1.5;

    printf("a = %0.4f\n", a);
    print_vector(4, x, "     x = ");
    print_vector(4, y, "     y = ");
    cblas_daxpy(4, a, x, 1, y, 1);
    print_vector(4, y, "ax + y = ");
}
