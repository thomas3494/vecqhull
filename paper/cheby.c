/* Adepted from Numerical Recipes, for analysing the Disk benchmark. */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double area(double x)
{
    return 2 * sin(x / 2) * (1 - cos(x / 2)) / (x - sin(x));
}

double *chebyshev_coef(double a, double b, double (*f)(double), int n)
{
    double *c = malloc(n * sizeof(double));

    double bma = 0.5 * (b - a);
    double bpa = 0.5 * (b + a);

    for (int j = 0; j < n; j++) {
        c[j] = 0.0;
        for (int k = 0; k < n; k++) {
            c[j] += f(cos(M_PI * (k + 0.5) / n) * bma + bpa) * 
                            cos(M_PI * j * (k + 0.5) / n);
        }
        c[j] *= 2.0 / n;
    }

    return c;
}

int main(int argc, char **argv)
{
    if (argc != 4) {
        printf("Usage: %s <DEGREE> <a> <b>\n", argv[0]);
        return EXIT_FAILURE;
    }

    int    n = atoi(argv[1]);
    double a = atof(argv[2]);
    double b = atof(argv[3]);

    double *c = chebyshev_coef(a, b, area, n);
    
    for (int i = 0; i < n; i++) {
        printf("c%d = %.17lf\n", i, c[i]);
    }

    free(c);
    return EXIT_SUCCESS;
}
