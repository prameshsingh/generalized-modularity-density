#define Inf_Sma 1e-18

struct graph
{
    int N;
    double n;
    int com;
    int **ml;
    int **nl;
    double **wl;
    double *dl;
};

struct part
{
    double Q;
    int com;
    int *pa;
};
