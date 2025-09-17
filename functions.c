#include "symnmf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* Implementation of the functions in the header file */

/* Calculate squared eucludean distance between two points */
double squared_euclidean_distance(double p[], double q[], int d)
{
    int i;
    double distance = 0.0;
    for (i = 0; i < d; ++i)
    {
        distance += pow((p[i] - q[i]), 2);
    }
    return distance;
}

/* Return the transpose of given matrix */
double **transpose(double **mat, int rows, int cols)
{
    int i, j;
    double **mat_t;
    mat_t = calloc(cols, sizeof(double *));
    if (mat_t == NULL)
    {
        return NULL;
    }
    for (i = 0; i < cols; i++)
    {
        mat_t[i] = calloc(rows, sizeof(double));
        if (mat_t[i] == NULL)
        {
            free(mat_t);
            return NULL;
        }
    }
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            mat_t[j][i] = mat[i][j];
        }
    }
    return mat_t;
}

/* return mat1*mat2 */
double **matrix_multiplication(double **mat1, double **mat2, int rows1, int rows2, int cols2)
{
    int i, j, k;
    double **result;
    if (mat1 == NULL || mat2 == NULL)
    {
        return NULL;
    }
    result = calloc(rows1, sizeof(double *));
    if (result == NULL)
    {
        return NULL;
    }
    for (i = 0; i < rows1; i++)
    {
        result[i] = calloc(cols2, sizeof(double));
        if (result[i] == NULL)
        {
            free(result);
            return NULL;
        }
    }
    /* Perform matrix multiplication */
    for (i = 0; i < rows1; i++)
    {
        for (k = 0; k < rows2; k++)
        {
            for (j = 0; j < cols2; j++)
            {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return result;
}

/* Calculate similarity matrix */
double **symC(double **data, int N, int d)
{
    int i, j;
    double **A;
    A = calloc(N, sizeof(double *));
    if (A == NULL)
    {
        return NULL;
    }
    for (i = 0; i < N; i++)
    {
        A[i] = calloc(N, sizeof(double));
        if (A[i] == NULL)
        {
            free(A);
            return NULL;
        }
    }
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (i != j)
            {
                A[i][j] = exp(-(squared_euclidean_distance(data[i], data[j], d) / 2));
            }
        }
    }
    return A;
}

/* Calculate diagonal degree matrix */
double ***ddgC(double **data, int N, int d)
{
    int i, j;
    double **D;
    double **A;
    double ***res;
    D = calloc(N, sizeof(double *));
    if (D == NULL)
    {
        return NULL;
    }
    for (i = 0; i < N; i++)
    {
        D[i] = calloc(N, sizeof(double));
        if (D[i] == NULL)
        {
            for (j = 0; j < i; j++)
                free(D[j]);
            free(D);
            return NULL;
        }
    }
    A = symC(data, N, d);
    if (A == NULL)
    {
        free(D);
        return NULL;
    }
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            D[i][i] += A[i][j];
        }
    }
    res = calloc(2, sizeof(double **));
    if (res == NULL)
    {
        for (i = 0; i < N; i++)
        {
            free(D[i]);
        }
        free(D);
        free(A);
        return NULL;
    }
    res[0] = A;
    res[1] = D;
    return res;
}

/* Calculate normalized similarity matrix */
double **normC(double **data, int N, int d)
{
    int i;
    double ***res = ddgC(data, N, d);
    double **A = res[0];
    double **D = res[1];
    double **W, **temp;
    if (A == NULL || D == NULL)
    {
        return NULL;
    }
    /* D = D^-0.5 */
    for (i = 0; i < N; i++)
    {
        if (D[i][i] != 0)
        {
            D[i][i] = pow(D[i][i], -0.5);
        }
    }
    /* W = D^-0.5 * A */
    W = matrix_multiplication(D, A, N, N, N);
    /* W = W * D^-0.5 = D^-0.5 * A * D^-0.5 */
    temp = matrix_multiplication(W, D, N, N, N);
    for (i = 0; i < N; i++)
    {
        free(W[i]);
        free(A[i]);
        free(D[i]);
    }
    free(A);
    free(D);
    free(res);
    free(W);
    W = temp;
    return W;
}

/* Calculate Frobenius norm */
double frobenius_norm(double **mat1, double **mat2, int rows, int col)
{
    double norm = 0.0;
    int i, j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < col; j++)
        {
            norm += pow(mat1[i][j] - mat2[i][j], 2);
        }
    }
    return norm;
}

/* Full symnmf algorithm */
double **symnmfC(double **H, double **W, int n, int k)
{
    double beta = 0.5;
    int i, j, iter = 0;
    double **newH, **numerator, **denominator, **H_t, **temp;
    newH = calloc(n, sizeof(double *));
    if (newH == NULL)
    {
        return NULL;
    }
    for (i = 0; i < n; i++)
    {
        newH[i] = calloc(k, sizeof(double));
        if (newH[i] == NULL)
        {
            free(newH);
            return NULL;
        }
    }
    /* Calculte H(1) */
    numerator = matrix_multiplication(W, H, n, n, k);
    H_t = transpose(H, n, k);
    if (H_t == NULL || numerator == NULL)
    {
        return NULL;
    }
    denominator = matrix_multiplication(H, H_t, n, k, n);
    denominator = matrix_multiplication(denominator, H, n, n, k);
    if (denominator == NULL)
    {
        return NULL;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < k; j++)
        {
            newH[i][j] = H[i][j] * (1 - beta + beta * (numerator[i][j] / denominator[i][j]));
        }
    }
    iter++;
    /* Iteratively calculate H(t) */
    while (frobenius_norm(newH, H, n, k) >= pow(10, -4) && iter < 300)
    {
        temp = H;
        H = newH;
        newH = temp;
        numerator = matrix_multiplication(W, H, n, n, k);
        H_t = transpose(H, n, k);
        if (H_t == NULL || numerator == NULL)
        {
            return NULL;
        }
        denominator = matrix_multiplication(H, H_t, n, k, n);
        denominator = matrix_multiplication(denominator, H, n, n, k);
        if (denominator == NULL)
        {
            return NULL;
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < k; j++)
            {
                newH[i][j] = H[i][j] * (1 - beta + beta * (numerator[i][j] / denominator[i][j]));
            }
        }
        iter++;
    }
    return newH;
}