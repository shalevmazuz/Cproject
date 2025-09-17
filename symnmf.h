#ifndef SYMNMF_H
#define SYMNMF_H

double **symC(double **data, int N, int d);
double ***ddgC(double **data, int N, int d);
double **normC(double **data, int N, int d);
double **symnmfC(double **H, double **W, int n, int k);

#endif