#include "symnmf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define COMMA 44

/* Calculate numbers of points */
int calc_N(char *file)
{
    int N = 0;
    char c;
    FILE *f;
    f = fopen(file, "r");
    if (f == NULL)
    {
        return -1;
    }
    for (c = getc(f); c != EOF; c = getc(f))
    {
        if (c == '\n') /* Increment count for every newline */
        {
            N++;
        }
    }
    fclose(f);
    return N;
}

/* Calculate dimension of points */
int calc_dim(char *file)
{
    int d = 1;
    char c;
    FILE *f;
    f = fopen(file, "r");
    if (f == NULL)
    {
        return -1;
    }
    for (c = getc(f); c != '\n'; c = getc(f))
    {
        if (c == COMMA) /* Increment count for every comma */
        {
            d++;
        }
    }
    fclose(f);
    return d;
}

int main(int argc, char *argv[])
{
    double **data;
    double element[1];
    double **sol;
    double ***sol2;
    int N, d, i, j, row = 0, col = 0, flag = 0;
    char *file = argv[2];
    char *goal = argv[1];
    FILE *f;
    if (argc < 1)
    {
        printf("An Error Has Occurred\n");
        return 1;
    }
    N = calc_N(file);
    d = calc_dim(file);
    if (N < 0 || d < 0) /* if calc_N or calc_dim failed */
    {
        printf("An Error Has Occurred\n");
        return 1;
    }
    f = fopen(file, "r");
    if (f == NULL)
    {
        printf("An Error Has Occurred\n");
        return 1;
    }
    data = calloc(N, sizeof(double *));
    if (data == NULL)
    {
        printf("An Error Has Occurred\n");
        return 1;
    }
    for (i = 0; i < N; i++)
    {
        data[i] = calloc(d, sizeof(double));
        if (data[i] == NULL)
        {
            printf("An Error Has Occurred\n");
            free(data);
            return 1;
        }
    }
    /* Read the file into 2d array */
    for (i = 0; fscanf(f, "%lf,", element) == 1; i++)
    {
        data[row][col] = element[0];
        if (col == d - 1)
        {
            row++;
            col = 0;
        }
        else
        {
            col++;
        }
    }
    fclose(f);
    if (strcmp(goal, "sym") == 0)
    {
        sol = symC(data, N, d);
        if (sol == NULL) /* if sym failed */
        {
            flag = 1;
        }
    }
    else if (strcmp(goal, "ddg") == 0)
    {
        sol2 = ddgC(data, N, d);
        if (sol2 == NULL) /* if ddg failed */
        {
            flag = 1;
        }
        sol = sol2[1];
        for (i = 0; i < N; i++)
        {
            free(sol2[0][i]);
        }
        free(sol2[0]);
        free(sol2);
    }
    else
    {
        sol = normC(data, N, d);
        if (sol == NULL) /* if norm failed */
        {
            flag = 1;
        }
    }
    if (flag == 1)
    {
        for (i = 0; i < N; i++)
        {
            free(data[i]);
        }
        free(data);
        printf("An Error Has Occurred\n");
        return 1;
    }
    /* Print required matrix */
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            printf("%.4f", sol[i][j]);
            if (j != N - 1)
            {
                printf(",");
            }
        }
        printf("\n");
    }
    for (i = 0; i < N; i++)
    {
        free(data[i]);
        free(sol[i]);
    }
    free(data);
    free(sol);
    return 0;
}