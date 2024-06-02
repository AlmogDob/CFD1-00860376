/* In this code every function and variable that have _f or _d in
it name it means that it uses float and double respectivly */

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define MAXDIR 100
#define MAXWORD 100
#define DEFAULT_N 100
#define PI 3.14159265359
#define x_d(i) start + (i) * h_d /* a macro for xi */
#define x_f(i) start + (i) * h_f

#define dprintSTRING(expr) printf(#expr " = %s\n", expr)    /* macro for easy debuging*/
#define dprintINT(expr) printf(#expr " = %d\n", expr)
#define dprintF(expr) printf(#expr " = %g\n", expr)
#define dprintD(expr) printf(#expr " = %g\n", expr)

void read_input(char *dir);
void allocat_and_fill_memmory(void);
void output_solution_d(char *dir, double *data);
void output_solution_f(char *dir, float *data);

float a_f(float xi);
float b_f(float xi);
float c_f(float xi);
float d_f(float xi);
double a_d(double xi);
double b_d(double xi);
double c_d(double xi);
double d_d(double xi);
float A_f(float xi);
float B_f(float xi);
float C_f(float xi);
float D_f(float xi);
double A_d(double xi);
double B_d(double xi);
double C_d(double xi);
double D_d(double xi);

void fill_matrix_d(double *A_vec_d, double *B_vec_d, double *C_vec_d, double *D_vec_d);
void fill_matrix_f(float *A_vec_f, float *B_vec_f, float *C_vec_f, float *D_vec_f);

int tridiag_f(float *a, float *b, float *c, float *d, float *u, int is, int ie);
int tridiag_d(double *a, double *b, double *c, double *d, double *u, int is, int ie);

double *A_vec_d = NULL, *B_vec_d = NULL, *C_vec_d = NULL, *D_vec_d = NULL, *u_vec_d = NULL; /* vectors for the system of equtions; double version*/
float *A_vec_f = NULL, *B_vec_f = NULL, *C_vec_f = NULL, *D_vec_f = NULL, *u_vec_f = NULL; /* vectors for the system of equtions; float version*/
int start = 0, end = 0, to_double = 1, N = -1,  /* interval: [start, end] */
to_dirichlet = 0, to_neumann = 0;   /* 'flags' for the typ of boundary conditions */
float h_f = 0.0f, Y_start_f = 0.0f, Y_end_f = 0.0f, /* Y_start = Y0; Y_end = YN */
Y_prim_start_f = 0.0f, Y_prim_end_f = 0.0f;
double h_d = 0.0, Y_start_d = 0.0, Y_end_d = 0.0,
Y_prim_start_d = 0.0, Y_prim_end_d = 0.0;

int main(int argc, char const *argv[])
{
    int success;

    /* Geting the input and output directories */
    char input_dir[MAXDIR], output_dir[MAXDIR];

    if (--argc != 2) {
        fprintf(stderr, "ERROR: not right usage\nUsage: main 'input dir' 'output dir'\n");
        return -1;
    }

    strncpy(input_dir, (*(++argv)), MAXDIR);
    strncpy(output_dir, (*(++argv)), MAXDIR);

    if (input_dir[MAXDIR-1] != '\0') {
        fprintf(stderr, "Error: input too long\n");
        return -1;
    }
    if (output_dir[MAXDIR-1] != '\0') {
        fprintf(stderr, "Error: output too long\n");
        return -1;
    }
    /*------------------------------------------------------------*/

    read_input(input_dir);

    // dprintINT(start);   /* checking if I read the input right */
    // dprintINT(end);
    // dprintINT(N);
    // dprintINT(to_double);
    // dprintD(h_d);
    // dprintINT(to_dirichlet);
    // dprintINT(to_neumann);
    // dprintD(Y_start_d);
    // dprintD(Y_end_d);
    // dprintD(Y_prim_start_d);
    // dprintD(Y_prim_end_d);
    /*------------------------------------------------------------*/

    allocat_and_fill_memmory();

    /*------------------------------------------------------------*/

    if (to_double){
        fill_matrix_d(A_vec_d, B_vec_d, C_vec_d, D_vec_d);
        if (to_dirichlet) {
            success = tridiag_d(A_vec_d, B_vec_d, C_vec_d, D_vec_d, u_vec_d, 1, N-1);
        } else if (to_neumann) {
            success = tridiag_d(A_vec_d, B_vec_d, C_vec_d, D_vec_d, u_vec_d, 0, N);
        }
        if (success == 0) {
            output_solution_d(output_dir, u_vec_d);
        } else {
            fprintf(stderr, "Error: something went wrong in the solver\n");
            return 1;
        }
    } else {
        fill_matrix_f(A_vec_f, B_vec_f, C_vec_f, D_vec_f);
        if (to_dirichlet) {
            success = tridiag_f(A_vec_f, B_vec_f, C_vec_f, D_vec_f, u_vec_f, 1, N-1);
        } else if (to_neumann) {
            success = tridiag_f(A_vec_f, B_vec_f, C_vec_f, D_vec_f, u_vec_f, 0, N);
        }
        if (success == 0) {
            output_solution_f(output_dir, u_vec_f);
        } else {
            fprintf(stderr, "Error: something went wrong in the solver\n");
            return 1;
        }

    }

    return 0;
}

/* sets 'flags' and variables according to the input file */
void read_input(char *dir)
{
    FILE *fp = fopen(dir, "rt");
    char current_word[MAXWORD];
    float temp, temp1;

    if (!fp) {
        fprintf(stderr, "Error opening file: %s\n", strerror(errno));        
    }

    while(fscanf(fp, "%s", current_word) != EOF) {
        if (!strcmp(current_word, "interval:")) {
            fscanf(fp, "%d %d", &start, &end);
        }
        if (!strcmp(current_word, "double?")) {
            fscanf(fp, "%s", current_word);
            if (!strcmp(current_word, "yes")) {
                to_double = 1;
            }
            if (!strcmp(current_word, "no")) {
                to_double = 0;
            }
        }
        if (!strcmp(current_word, "N:")) {
            fscanf(fp, "%d ", &N);
        }
        if (!strcmp(current_word, "dirichlet?")) {
            fscanf(fp, "%s", current_word);
            if (!strcmp(current_word, "yes")) {
                to_dirichlet = 1;
            }
            if (!strcmp(current_word, "no")) {
                to_dirichlet = 0;
            }
            if (to_double) {
                fscanf(fp, "%g %g", &temp, &temp1);
                Y_start_d = (double)temp;
                Y_end_d = (double)temp1;
            } else {
                fscanf(fp, "%e %e", &Y_start_f, &Y_end_f);
            }
        }
        if (!strcmp(current_word, "neumann?")) {
            fscanf(fp, "%s", current_word);
            if (!strcmp(current_word, "yes")) {
                to_neumann = 1;
            }
            if (!strcmp(current_word, "no")) {
                to_neumann = 0;
            }
            if (to_double) {
                fscanf(fp, "%g %g", &temp, &temp1);
                Y_prim_start_d = (double)temp;
                Y_prim_end_d = (double)temp1;
            } else {
                fscanf(fp, "%e %e", &Y_prim_start_f, &Y_prim_end_f);
            }
        }
    }
    if (!to_neumann && !to_dirichlet) {
        fprintf(stderr, "Error: unknow boundary conditions\n");
        exit(1);
    }
    if (N == -1) {
        N = DEFAULT_N;
    }

    /* Calculation h*/
    if (to_double) {
        h_d = (end - start)/(double)N;
    } else {
        h_f = (end - start)/(float)N;
    }

}

/* allocating the vectors and filling with zero */
void allocat_and_fill_memmory(void)
{
    int i;
    A_vec_d = (double *)malloc(sizeof(double) * (N + 1));
    for (i = 0; i < N+1; i++) {
        A_vec_d[i] = 0;
    }
    B_vec_d = (double *)malloc(sizeof(double) * (N + 1));
    for (i = 0; i < N+1; i++) {
        B_vec_d[i] = 0;
    }
    C_vec_d = (double *)malloc(sizeof(double) * (N + 1));
    for (i = 0; i < N+1; i++) {
        C_vec_d[i] = 0;
    }
    D_vec_d = (double *)malloc(sizeof(double) * (N + 1));
    for (i = 0; i < N+1; i++) {
        D_vec_d[i] = 0;
    }
    u_vec_d = (double *)malloc(sizeof(double) * (N + 1));
    for (i = 0; i < N+1; i++) {
        u_vec_d[i] = 0;
    }
    A_vec_f = (float *)malloc(sizeof(float) * (N + 1));
    for (i = 0; i < N+1; i++) {
        A_vec_f[i] = 0;
    }
    B_vec_f = (float *)malloc(sizeof(float) * (N + 1));
    for (i = 0; i < N+1; i++) {
        B_vec_f[i] = 0;
    }
    C_vec_f = (float *)malloc(sizeof(float) * (N + 1));
    for (i = 0; i < N+1; i++) {
        C_vec_f[i] = 0;
    }
    D_vec_f = (float *)malloc(sizeof(float) * (N + 1));
    for (i = 0; i < N+1; i++) {
        D_vec_f[i] = 0;
    }
    u_vec_f = (float *)malloc(sizeof(float) * (N + 1));
    for (i = 0; i < N+1; i++) {
        u_vec_f[i] = 0;
    }
}

/* output data; double version */
void output_solution_d(char *dir, double *data)
{
    FILE *fp = fopen(dir, "wt");
    int i;
    char buff[20];

    /* with dirichlet boundary conditions the output has elements from 1 to N-1 */
    if (to_dirichlet) {
        for (i = 1; i < N; i++) {
            // sprintf(buff, "%%%dd", (int)log10(N));
            // fprintf(fp, buff, i);
            // fprintf(fp, ", %g\n", data[i]);
            fprintf(fp, "%g\n", data[i]);
        }
    }

    /* with neumann boundary conditions the output has elements from 0 to N */
    if (to_neumann) {
        for (i = 0; i < N+1; i++) {
            // sprintf(buff, "%%%dd", (int)log10(N));
            // fprintf(fp, buff, i);
            // fprintf(fp, ", %g\n", data[i]);
            fprintf(fp, "%g\n", data[i]);
        }
    }
    fclose(fp);
    printf("Output %d lines to %s in the interval: [%d, %d]\n", N, dir, start, end);
    printf("Start conditions: ");
    if (to_dirichlet) {
        printf("dirichlet Y0 = %g, YN = %g\n", Y_start_d, Y_end_d);
    } else {
        printf("neumann Y'0 = %g, Y'N = %g\n", Y_prim_start_d, Y_prim_end_d);
    }
    printf("Type: ");
    printf("double\n");

}

/* output data; float version */
void output_solution_f(char *dir, float *data)
{
    FILE *fp = fopen(dir, "wt");
    int i;
    char buff[20];

    /* with dirichlet boundary conditions the output has elements from 1 to N-1 */
    if (to_dirichlet) {
        for (i = 1; i < N; i++) {
            // sprintf(buff, "%%%dd", (int)log10(N));
            // fprintf(fp, buff, i);
            // fprintf(fp, ", %g\n", data[i]);
            fprintf(fp, "%g\n", data[i]);
        }
    }

    /* with neumann boundary conditions the output has elements from 0 to N */
    if (to_neumann) {
        for (i = 0; i < N+1; i++) {
            // sprintf(buff, "%%%dd", (int)log10(N));
            // fprintf(fp, buff, i);
            // fprintf(fp, ", %g\n", data[i]);
            fprintf(fp, "%g\n", data[i]);
        }
    }
    fclose(fp);
    printf("Output %d lines to %s in the interval: [%d, %d]\n", N, dir, start, end);
    printf("Start conditions: ");
    if (to_dirichlet) {
        printf("dirichlet Y0 = %g, YN = %g\n", Y_start_f, Y_end_f);
    } else {
        printf("neumann Y'0 = %g, Y'N = %g\n", Y_prim_start_f, Y_prim_end_f);
    }
    printf("Type: ");
    printf("float\n");

}

/* The function a in float*/
float a_f(float xi)
{
    xi = 1;
    return xi;
}

/* The function b in float*/
float b_f(float xi)
{
    return xi;
}

/* The function c in float*/
float c_f(float xi)
{
    return xi*xi;
}

/* The function d in float*/
float d_f(float xi)
{
    return sinf(2*PI*xi)+ cosf(2*PI*xi);    
}

/* The function a in double*/
double a_d(double xi)
{
    xi = 1;
    return xi;
}

/* The function b in double*/
double b_d(double xi)
{
    return xi;
}

/* The function c in double*/
double c_d(double xi)
{
    return xi*xi;
}

/* The function d in double*/
double d_d(double xi)
{
    // printf("d_d = %g\n", sin(2*PI*xi)+ cos(2*PI*xi));    
    return sin(2*PI*xi)+ cos(2*PI*xi);    
}

/* A from the discrete from of the equation in float */
float A_f(float xi)
{
    return (a_f(xi)/(h_f*h_f) - b_f(xi)/(2*h_f));
}

/* B from the discrete from of the equation in float */
float B_f(float xi)
{
    return (-2*a_f(xi)/(h_f*h_f) +c_f(xi));
}

/* C from the discrete from of the equation in float */
float C_f(float xi)
{
    return (a_f(xi)/(h_f*h_f) +b_f(xi)/(2*h_f));
}

/* D from the discrete from of the equation in float */
float D_f(float xi)
{
    return d_f(xi);
}

/* A from the discrete from of the equation in double */
double A_d(double xi)
{
    return (a_d(xi)/(h_d*h_d) - b_d(xi)/(2*h_d));
}

/* B from the discrete from of the equation in double */
double B_d(double xi)
{
    return (-2*a_d(xi)/(h_d*h_d) +c_d(xi));
}

/* C from the discrete from of the equation in double */
double C_d(double xi)
{
    return (a_d(xi)/(h_d*h_d) +b_d(xi)/(2*h_d));
}

/* D from the discrete from of the equation in double */
double D_d(double xi)
{
    return d_d(xi);
}

/* This function creats the vectors needed to call the solver; double version */
void fill_matrix_d(double *A_vec_d, double *B_vec_d, double *C_vec_d, double *D_vec_d)
{
    int i;

    if (to_dirichlet) {
        for (i = 1; i < N+1; i++) { /* going over all the points */
            if (i != N) {   /* for dirichlet we don't need the N's index */
                /* LHS */
                if (i != 1) {   /* there is no A1 in dirichlet */
                    A_vec_d[i] = A_d(x_d(i));
                }
                B_vec_d[i] = B_d(x_d(i));
                if (i != N-1) {
                    C_vec_d[i] = C_d(x_d(i));
                }
                /* RHS */
                D_vec_d[1] = D_d(x_d(1)) - A_d(x_d(1)) * Y_start_d;
                if (i > 1 && i < N-1) {
                    D_vec_d[i] = D_d(x_d(i));
                }
                D_vec_d[N-1] = D_d(x_d(N-1)) - C_d(x_d(N-1)) * Y_end_d;
            }                
        }
    }

    if (to_neumann) {
        for (i = 0; i < N+1; i++) { /* going over all the points */
            /* LHS */
            if (i != 0) {   /* there is no A_vec[0] in neumann */
                if (i == N) {   /* A_vec[N] is different */
                    A_vec_d[i] = A_d(x_d(i)) + C_d(x_d(i));
                } else {
                    A_vec_d[i] = A_d(x_d(i));
                }
            }
            B_vec_d[i] = B_d(x_d(i));
            if (i != N) {   /* there is no A_vec[N] in neumann */
                if (i == 0) {   /* C_vec[0] is different */
                    C_vec_d[i] = A_d(x_d(i)) + C_d(x_d(i));
                } else {
                    C_vec_d[i] = C_d(x_d(i));
                }
            }
            /* RHS */
            D_vec_d[0] = D_d(x_d(0)) + 2 * h_d * A_d(x_d(0)) * Y_prim_start_d;
            if (i > 0 && i < N) {
                D_vec_d[i] = D_d(x_d(i));
            }
            D_vec_d[N] = D_d(x_d(N)) - 2 * h_d * C_d(x_d(N)) * Y_prim_end_d;
        }
    }
}

/* This function creats the vectors needed to call the solver; float version */
void fill_matrix_f(float *A_vec_f, float *B_vec_f, float *C_vec_f, float *D_vec_f)
{
    int i;

    if (to_dirichlet) {
        for (i = 1; i < N+1; i++) { /* going over all the points */
            if (i != N) {   /* for dirichlet we don't need the N's index */
                /* LHS */
                if (i != 1) {   /* there is no A1 in dirichlet */
                    A_vec_f[i] = A_f(x_f(i));
                }
                B_vec_f[i] = B_f(x_f(i));
                if (i != N-1) {
                    C_vec_f[i] = C_f(x_f(i));
                }
                /* RHS */
                D_vec_f[1] = D_f(x_f(1)) - A_f(x_f(1)) * Y_start_f;
                if (i > 1 && i < N-1) {
                    D_vec_f[i] = D_f(x_f(i));
                }
                D_vec_f[N-1] = D_f(x_f(N-1)) - C_f(x_f(N-1)) * Y_end_f;
            }                
        }
    }

    if (to_neumann) {
        for (i = 0; i < N+1; i++) { /* going over all the points */
            /* LHS */
            if (i != 0) {   /* there is no A_vec[0] in neumann */
                if (i == N) {   /* A_vec[N] is different */
                    A_vec_f[i] = A_f(x_f(i)) + C_f(x_f(i));
                } else {
                    A_vec_f[i] = A_f(x_f(i));
                }
            }
            B_vec_f[i] = B_f(x_f(i));
            if (i != N) {   /* there is no A_vec[N] in neumann */
                if (i == 0) {   /* C_vec[0] is different */
                    C_vec_f[i] = A_f(x_f(i)) + C_f(x_f(i));
                } else {
                    C_vec_f[i] = C_f(x_f(i));
                }
            }
            /* RHS */
            D_vec_f[0] = D_f(x_f(0)) + 2 * h_f * A_f(x_f(i)) * Y_prim_start_f;
            if (i > 0 && i < N) {
                D_vec_f[i] = D_f(x_f(i));
            }
            D_vec_f[N] = D_f(x_f(N)) - 2 * h_f * C_f(x_f(N)) * Y_prim_end_f;
        }
    }
}

/* a, b, c, are the vectors of the diagonal and the two
off-diagonals. The vector d is the RHS vector, the vector
u is the solution vector, "is" is the starting point, and
ie is the last point; float version */
int tridiag_f(float *a, float *b, float *c, float *d, float *u, int is, int ie)
{
    int i;
    float beta;

    /*test*/
    dprintF(b[ie]);
    dprintF(c[ie-1]);
    dprintF(b[ie] - c[ie - 1] * (-1.999900));
    /*test*/
    for (i = is + 1; i <= ie; i++)
    {
        if (b[i - 1] == 0.)
            return (1);
        beta = a[i] / b[i - 1];
        /*test*/
        if (i == ie) {
            dprintF(b[i]);
            dprintF(c[i-1]);
            dprintF(beta);
            dprintF(b[i] - c[i - 1] * beta);
            dprintF(d[i]);
        }
        /*test*/
        b[i] = b[i] - c[i - 1] * beta;
        d[i] = d[i] - d[i - 1] * beta;
        
    }

    u[ie] = d[ie] / b[ie];
    for (i = ie - 1; i >= is; i--)
    {
        u[i] = (d[i] - c[i] * u[i + 1]) / b[i];
    }
    return (0);
}

/* a, b, c, are the vectors of the diagonal and the two
off-diagonals. The vector d is the RHS vector, the vector
u is the solution vector, "is" is the starting point, and
ie is the last point; double version */
int tridiag_d(double *a, double *b, double *c, double *d, double *u, int is, int ie)
{
    int i;
    double beta;

    for (i = is + 1; i <= ie; i++)
    {
        if (b[i - 1] == 0.)
            return (1);
        beta = a[i] / b[i - 1];
        b[i] = b[i] - c[i - 1] * beta;
        d[i] = d[i] - d[i - 1] * beta;
    }

    u[ie] = d[ie] / b[ie];
    for (i = ie - 1; i >= is; i--)
    {
        u[i] = (d[i] - c[i] * u[i + 1]) / b[i];
    }
    return (0);
}