#include <stdio.h>

int main(int argc, char **argv)
{
    if (argc != 1) {
        printf("Usage: %s\n"
               "\tTakes points in double precision on stdin, returns"
               "\tthem in float to stdout\n", argv[0]);
        return EXIT_FAILURE;
    }

    double x;
    float xf;
    while (fread(&x, sizeof(double), 1, stdin) == 1) {
        xf = (float)x;
        fwrite(&xf, sizeof(float), 1, stdout);
    }

    return EXIT_SUCCESS;
}
