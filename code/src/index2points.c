#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef struct {
    double x;
    double y;
} Point;

Point *parse_point2d(size_t *n, FILE *points_f)
{
    size_t line_length;
    char *file_header = NULL;

    getline(&file_header, &line_length, points_f);

    if (strcmp(file_header, "pbbs_sequencePoint2d\n") != 0) {
        printf("Error: not in pbbs_sequencePoint2d format\n");
        *n = 0;
        return NULL;
    }

    free(file_header);

    size_t size = 1;
    Point *P = (Point *)malloc(size * sizeof(Point));

    size_t i = 0;
    while (fscanf(points_f, "%lf %lf", &P[i].x, &P[i].y) == 2) {
        i++;
        if (i >= size) {
            size *= 2;
            P = (Point *)realloc(P, size * sizeof(Point));
        }
    }

    *n = i;
    P = (Point *)realloc(P, i * sizeof(Point));

    return P;
}

size_t *parse_sequenceInt(size_t *n, FILE *indices_f)
{
    size_t line_length;
    char *file_header = NULL;

    getline(&file_header, &line_length, indices_f);

    if (strcmp(file_header, "sequenceInt\n") != 0) {
        printf("Error: not in sequenceInt format\n");
        *n = 0;
        return NULL;
    }

    free(file_header);

    size_t size = 1;
    size_t *indices = (size_t *)malloc(size * sizeof(size_t));

    size_t i = 0;
    while (fscanf(indices_f, "%zu", &indices[i]) == 1) {
        i++;
        if (i >= size) {
            size *= 2;
            indices = (size_t *)realloc(indices, size * sizeof(size_t));
        }
    }

    *n = i;
    indices = (size_t *)realloc(indices, i * sizeof(size_t));

    return indices;
}

int main(int argc, char **argv)
{
    if (argc != 3) {
        printf("Usage: %s POINTS INDICES\n"
               "\tPOINTS: full set if points in pbbs_sequencePoint2d format.\n"
               "\tINDICES: full set if indices in sequenceInt format.\n"
               "\tOutputs the subset of points to stdout"
               " in pbbs_sequencePoint2d format\n", argv[0]);
        return EXIT_FAILURE;
    }

    size_t n_p;
    FILE *points_f = fopen(argv[1], "r");
    if (points_f == NULL) {
        printf("%s needs to be a readable file\n", argv[1]);
        return EXIT_FAILURE;
    }
    Point *P = parse_point2d(&n_p, points_f);
    fclose(points_f);
    if (P == NULL) {
        printf("Error reading points\n");
        return EXIT_FAILURE;
    }

    size_t n_i;
    FILE *indices_f = fopen(argv[2], "r");
    if (indices_f == NULL) {
        free(P);
        printf("%s needs to be a readable file\n", argv[2]);
        abort();
    }
    size_t *indices = parse_sequenceInt(&n_i, indices_f);
    fclose(indices_f);
    if (indices == NULL) {
        printf("Erorr reading indices\n");
        return EXIT_FAILURE;
    }

    printf("pbbs_sequencePoint2d\n");
    for (size_t i = 0; i < n_i; i++) {
        size_t j = indices[i];
        printf("%.17e %.17e\n", P[j].x, P[j].y);
    }

    free(P);
    free(indices);
}
