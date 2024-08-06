#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

int main(int argc, char **argv)
{
    if (argc != 1) {
        printf("Usage: %s\n"
                "\tTakes points on stdin in text format returns in binary\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    size_t line_length;
    char *file_header = NULL;
    if (getline(&file_header, &line_length, stdin) == -1) {
        perror("Reading line failed\n");
    }

    if (strcmp(file_header, "pbbs_sequencePoint2d\n") != 0) {
        printf("Error: not in pbbs_sequencePoint2d format\n");
        abort();
    }

    double x, y;
    while (scanf("%lf %lf\n", &x, &y) == 2) {
        fwrite(&x, sizeof(double), 1, stdout);
        fwrite(&y, sizeof(double), 1, stdout);
    }

    free(file_header);

    return EXIT_SUCCESS;
}
