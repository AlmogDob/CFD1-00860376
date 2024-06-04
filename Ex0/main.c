/* Example for input file (The order does not matter):
```
interval:
0 1
double?
yes
N:
100
dirichlet?
yes
0.0 1.0
neumann?
no
1.0 -1.0
output_check_solution?
yes
```
Program call: ./main 'input directory' 'output directory' */

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
void boundary_conditions_d(double *A_vec_d, double *B_vec_d, double *C_vec_d, double *D_vec_d);
void boundary_conditions_f(float *A_vec_f, float *B_vec_f, float *C_vec_f, float *D_vec_f);

int tridiag_f(float *a, float *b, float *c, float *d, float *u, int is, int ie);
int tridiag_d(double *a, double *b, double *c, double *d, double *u, int is, int ie);

void check_solution_d(double *check_diff_d, double *A_vec_d, double *B_vec_d, double *C_vec_d, double *D_vec_d, double *u_vec_d);
void check_solution_f(float *check_diff_f, float *A_vec_f, float *B_vec_f, float *C_vec_f, float *D_vec_f, float *u_vec_f);

int start = 0, end = 0, to_double = 1, N = -1,  /* interval: [start, end] */
to_dirichlet = 0, to_neumann = 0, to_check = 0;   /* 'flags' for the typ of boundary conditions and checking the solution */
float h_f = 0.0f, Y_start_f = 0.0f, Y_end_f = 0.0f, /* Y_start = Y0; Y_end = YN */
Y_prim_start_f = 0.0f, Y_prim_end_f = 0.0f;
double h_d = 0.0, Y_start_d = 0.0, Y_end_d = 0.0,
Y_prim_start_d = 0.0, Y_prim_end_d = 0.0;

int main(int argc, char const *argv[])
{
    int success, i;

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
    // dprintINT(to_check);
    // dprintD(Y_start_d);
    // dprintD(Y_end_d);
    // dprintD(Y_prim_start_d);
    // dprintD(Y_prim_end_d);
    /*------------------------------------------------------------*/

    /* vector memmory allocation and filling with zero */
    double *A_vec_d = (double *)malloc(sizeof(double) * (N + 1));
    for (i = 0; i < N+1; i++) {
        A_vec_d[i] = 0;
    }
    double *B_vec_d = (double *)malloc(sizeof(double) * (N + 1));
    for (i = 0; i < N+1; i++) {
        B_vec_d[i] = 0;
    }
    double *C_vec_d = (double *)malloc(sizeof(double) * (N + 1));
    for (i = 0; i < N+1; i++) {
        C_vec_d[i] = 0;
    }
    double *D_vec_d = (double *)malloc(sizeof(double) * (N + 1));
    for (i = 0; i < N+1; i++) {
        D_vec_d[i] = 0;
    }
    double *u_vec_d = (double *)malloc(sizeof(double) * (N + 1));
    for (i = 0; i < N+1; i++) {
        u_vec_d[i] = 0;
    }
    double *check_diff_d = (double *)malloc(sizeof(double) * (N + 1));
    for (i = 0; i < N+1; i++) {
        check_diff_d[i] = 0;
    }
    float *A_vec_f = (float *)malloc(sizeof(float) * (N + 1));
    for (i = 0; i < N+1; i++) {
        A_vec_f[i] = 0;
    }
    float *B_vec_f = (float *)malloc(sizeof(float) * (N + 1));
    for (i = 0; i < N+1; i++) {
        B_vec_f[i] = 0;
    }
    float *C_vec_f = (float *)malloc(sizeof(float) * (N + 1));
    for (i = 0; i < N+1; i++) {
        C_vec_f[i] = 0;
    }
    float *D_vec_f = (float *)malloc(sizeof(float) * (N + 1));
    for (i = 0; i < N+1; i++) {
        D_vec_f[i] = 0;
    }
    float *u_vec_f = (float *)malloc(sizeof(float) * (N + 1));
    for (i = 0; i < N+1; i++) {
        u_vec_f[i] = 0;
    }
    float *check_diff_f = (float *)malloc(sizeof(float) * (N + 1));
    for (i = 0; i < N+1; i++) {
        check_diff_f[i] = 0;
    }

    /*------------------------------------------------------------*/

    if (to_double){
        fill_matrix_d(A_vec_d, B_vec_d, C_vec_d, D_vec_d);
        if (to_dirichlet) {
            success = tridiag_d(A_vec_d, B_vec_d, C_vec_d, D_vec_d, u_vec_d, 1, N-1);
        } else if (to_neumann) {
            success = tridiag_d(A_vec_d, B_vec_d, C_vec_d, D_vec_d, u_vec_d, 0, N);
        }
        if (success == 0) {
            if (!to_check) {     /* if we are not chekc the solution */
                output_solution_d(output_dir, u_vec_d);
            }
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
            if (!to_check) {     /* if we are not chekc the solution */
                output_solution_f(output_dir, u_vec_f);
            }
        } else {
            fprintf(stderr, "Error: something went wrong in the solver\n");
            return 1;
        }
    }

    if (to_double) {
        check_solution_d(check_diff_d, A_vec_d, B_vec_d, C_vec_d, D_vec_d, u_vec_d);
        if (to_check) {
            output_solution_d(output_dir, check_diff_d);
        }
    } else {
        check_solution_f(check_diff_f, A_vec_f, B_vec_f, C_vec_f, D_vec_f, u_vec_f);
        if (to_check) {
            output_solution_f(output_dir, check_diff_f);
        }
    }
    

    return 0;
}

/* sets 'flags' and variables according to the input file
argument list: dir - the directory of the input file */
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
        if (!strcmp(current_word, "output_check_solution?")) {
            fscanf(fp, "%s", current_word);
            if (!strcmp(current_word, "yes")) {
                to_check = 1;
            }
            if (!strcmp(current_word, "no")) {
                to_check = 0;
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

/* output data; double version
argument list: dir - the directory of the output file.
data - the solution vector */
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

/* output data; float version
argument list: dir - the directory of the output file.
data - the solution vector */
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

/* The function a in float
argument list: xi - the point where we want the value of the function */
float a_f(float xi)
{
    xi = 1;
    return xi;
}

/* The function b in float
argument list: xi - the point where we want the value of the function */
float b_f(float xi)
{
    return xi;
}

/* The function c in float
argument list: xi - the point where we want the value of the function */
float c_f(float xi)
{
    return xi*xi;
}

/* The function d in float
argument list: xi - the point where we want the value of the function */
float d_f(float xi)
{
    return sinf(2*PI*xi)+ cosf(2*PI*xi);    
}

/* The function a in double
argument list: xi - the point where we want the value of the function */
double a_d(double xi)
{
    xi = 1;
    return xi;
}

/* The function b in double
argument list: xi - the point where we want the value of the function */
double b_d(double xi)
{
    return xi;
}

/* The function c in double
argument list: xi - the point where we want the value of the function */
double c_d(double xi)
{
    return xi*xi;
}

/* The function d in double
argument list: xi - the point where we want the value of the function */
double d_d(double xi)
{
    // printf("d_d = %g\n", sin(2*PI*xi)+ cos(2*PI*xi));    
    return sin(2*PI*xi)+ cos(2*PI*xi);    
}

/* A from the discrete form of the equation in float
argument list: xi - the point where we want the value of the function */
float A_f(float xi)
{
    return (a_f(xi)/(h_f*h_f) - b_f(xi)/(2*h_f));
}

/* B from the discrete form of the equation in float
argument list: xi - the point where we want the value of the function */
float B_f(float xi)
{
    return (-2*a_f(xi)/(h_f*h_f) +c_f(xi));
}

/* C from the discrete form of the equation in float
argument list: xi - the point where we want the value of the function */
float C_f(float xi)
{
    return (a_f(xi)/(h_f*h_f) +b_f(xi)/(2*h_f));
}

/* D from the discrete form of the equation in float
argument list: xi - the point where we want the value of the function */
float D_f(float xi)
{
    return d_f(xi);
}

/* A from the discrete form of the equation in double
argument list: xi - the point where we want the value of the function */
double A_d(double xi)
{
    return (a_d(xi)/(h_d*h_d) - b_d(xi)/(2*h_d));
}

/* B from the discrete form of the equation in double
argument list: xi - the point where we want the value of the function */
double B_d(double xi)
{
    return (-2*a_d(xi)/(h_d*h_d) +c_d(xi));
}

/* C from the discrete form of the equation in double
argument list: xi - the point where we want the value of the function */
double C_d(double xi)
{
    return (a_d(xi)/(h_d*h_d) +b_d(xi)/(2*h_d));
}

/* D from the discrete form of the equation in double
argument list: xi - the point where we want the value of the function */
double D_d(double xi)
{
    return d_d(xi);
}

/* This function creats the vectors needed to call the solver; double version
argument list: A_vec_d, B_vec_d, C_vec_d - are the diagonals of the LHS matrix.
               D_vec_d - the RHS */
void fill_matrix_d(double *A_vec_d, double *B_vec_d, double *C_vec_d, double *D_vec_d)
{
    int i;

    if (to_dirichlet) {
        for (i = 1; i < N; i++) { /* going over all the points for dirichlet */
            /* LHS */
            if (i != 1) {   /* there is no A1 in dirichlet */
                A_vec_d[i] = A_d(x_d(i));
            }
            B_vec_d[i] = B_d(x_d(i));
            if (i != N-1) { /* there is no C_N-1 in dirichlet*/
                C_vec_d[i] = C_d(x_d(i));
            }
            /* RHS */
            D_vec_d[i] = D_d(x_d(i));
            boundary_conditions_d(A_vec_d, B_vec_d, C_vec_d, D_vec_d);
        }
    }

    if (to_neumann) {
        for (i = 0; i < N+1; i++) { /* going over all the points */
            /* LHS */
            if (i != 0) {   /* there is no A_vec[0] in neumann */
                A_vec_d[i] = A_d(x_d(i));
            }
            B_vec_d[i] = B_d(x_d(i));
            if (i != N) {   /* there is no C_vec[N] in neumann */
                C_vec_d[i] = C_d(x_d(i));
            }
            /* RHS */
            D_vec_d[i] = D_d(x_d(i));
            boundary_conditions_d(A_vec_d, B_vec_d, C_vec_d, D_vec_d);
        }
    }
}

/* This function creats the vectors needed to call the solver; float version
argument list: A_vec_d, B_vec_d, C_vec_d - are the diagonals of the LHS matrix.
               D_vec_d - the RHS */
void fill_matrix_f(float *A_vec_f, float *B_vec_f, float *C_vec_f, float *D_vec_f)
{
    int i;

    if (to_dirichlet) {
        for (i = 1; i < N; i++) { /* going over all the points for dirichlet */
            /* LHS */
            if (i != 1) {   /* there is no A1 in dirichlet */
                A_vec_f[i] = A_f(x_f(i));
            }
            B_vec_f[i] = B_f(x_f(i));
            if (i != N-1) { /* there is no C_N-1 in dirichlet*/
                C_vec_f[i] = C_f(x_f(i));
            }
            /* RHS */
            D_vec_f[i] = D_f(x_f(i));
            boundary_conditions_f(A_vec_f, B_vec_f, C_vec_f, D_vec_f);
        }
    }

    if (to_neumann) {
        for (i = 0; i < N+1; i++) { /* going over all the points */
            /* LHS */
            if (i != 0) {   /* there is no A_vec[0] in neumann */
                A_vec_f[i] = A_f(x_f(i));
            }
            B_vec_f[i] = B_f(x_f(i));
            if (i != N) {   /* there is no C_vec[N] in neumann */
                C_vec_f[i] = C_f(x_f(i));
            }
            /* RHS */
            D_vec_f[i] = D_f(x_f(i));
            boundary_conditions_f(A_vec_f, B_vec_f, C_vec_f, D_vec_f);
        }
    }
}

/* Set the boundary conditions in the vectors; double version
argument list: A_vec_d, B_vec_d, C_vec_d - are the diagonals of the LHS matrix.
               D_vec_d - the RHS */
void boundary_conditions_d(double *A_vec_d, double *B_vec_d, double *C_vec_d, double *D_vec_d)
{
    if (to_dirichlet) {
        D_vec_d[1] = D_d(x_d(1)) - A_d(x_d(1)) * Y_start_d;
        D_vec_d[N-1] = D_d(x_d(N-1)) - C_d(x_d(N-1)) * Y_end_d;
    }
    if (to_neumann) {
        A_vec_d[N] = A_d(x_d(N)) + C_d(x_d(N));
        C_vec_d[0] = A_d(x_d(0)) + C_d(x_d(0));
        D_vec_d[0] = D_d(x_d(0)) + 2 * h_d * A_d(x_d(0)) * Y_prim_start_d;
        D_vec_d[N] = D_d(x_d(N)) - 2 * h_d * C_d(x_d(N)) * Y_prim_end_d;
    }
    (void)B_vec_d;
}

/* Set the boundary conditions in the vectors; float version
argument list: A_vec_f, B_vec_f, C_vec_f - are the diagonals of the LHS matrix.
               D_vec_f - the RHS */
void boundary_conditions_f(float *A_vec_f, float *B_vec_f, float *C_vec_f, float *D_vec_f)
{
    if (to_dirichlet) {
        D_vec_f[1] = D_f(x_f(1)) - A_f(x_f(1)) * Y_start_f;
        D_vec_f[N-1] = D_f(x_f(N-1)) - C_f(x_f(N-1)) * Y_end_f;
    }
    if (to_neumann) {
        A_vec_f[N] = A_f(x_f(N)) + C_f(x_f(N));
        C_vec_f[0] = A_f(x_f(0)) + C_f(x_f(0));
        D_vec_f[0] = D_f(x_f(0)) + 2 * h_f * A_f(x_f(0)) * Y_prim_start_f;
        D_vec_f[N] = D_f(x_f(N)) - 2 * h_f * C_f(x_f(N)) * Y_prim_end_f;
    }
    (void)B_vec_f;
}

/* a, b, c, are the vectors of the diagonal and the two
off-diagonals. The vector d is the RHS vector, the vector
u is the solution vector, "is" is the starting point, and
ie is the last point; float version */
int tridiag_f(float *a, float *b, float *c, float *d, float *u, int is, int ie)
{
    int i;
    float beta;

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

/* Calcultas Ax-B and puts the answer in check_diff; double version
argument list: A_vec_d, B_vec_d, C_vec_d - LHS.
D_vec_d - RHS.
u_vec_d - the solution vector.
check_diff_d - the vector of the differences. */
void check_solution_d(double *check_diff_d, double *A_vec_d,
                      double *B_vec_d, double *C_vec_d,
                      double *D_vec_d, double *u_vec_d)
{
    int i;

    fill_matrix_d(A_vec_d, B_vec_d, C_vec_d,D_vec_d);
    if (to_dirichlet) {
        check_diff_d[1] = B_vec_d[1] * u_vec_d[1] + C_vec_d[1] * u_vec_d[2] - D_vec_d[1];
        for (i = 2; i < N-1; i++) {
            check_diff_d[i] = A_vec_d[i] * u_vec_d[i-1] + B_vec_d[i] * u_vec_d[i] + C_vec_d[i] * u_vec_d[i+1] - D_vec_d[i];
        }
        check_diff_d[N-1] = A_vec_d[N-1] * u_vec_d[N-2] + B_vec_d[N-1] * u_vec_d[N-1] - D_vec_d[N-1];
    }
    if (to_neumann) {
        check_diff_d[0] = B_vec_d[0] * u_vec_d[0] + (A_vec_d[0] + C_vec_d[0]) * u_vec_d[1] - D_vec_d[0];
        for (i = 1; i < N; i++) {
            check_diff_d[i] = A_vec_d[i] * u_vec_d[i-1] + B_vec_d[i] * u_vec_d[i] + C_vec_d[i] * u_vec_d[i+1] - D_vec_d[i];
        }
        check_diff_d[N] = (A_vec_d[N] + C_vec_d[N]) * u_vec_d[N-1] + B_vec_d[N] * u_vec_d[N] - D_vec_d[N];
    }
}

/* Calcultas Ax-B and puts the answer in check_diff; float version
argument list: A_vec_f, B_vec_f, C_vec_f - LHS.
D_vec_f - RHS.
u_vec_f - the solution vector.
check_diff_f - the vector of the differences. */
void check_solution_f(float *check_diff_f, float *A_vec_f,
                      float *B_vec_f, float *C_vec_f,
                      float *D_vec_f, float *u_vec_f)
{
    int i;

    fill_matrix_f(A_vec_f, B_vec_f, C_vec_f,D_vec_f);
    if (to_dirichlet) {
        check_diff_f[1] = B_vec_f[1] * u_vec_f[1] + C_vec_f[1] * u_vec_f[2] - D_vec_f[1];
        for (i = 2; i < N-1; i++) {
            check_diff_f[i] = A_vec_f[i] * u_vec_f[i-1] + B_vec_f[i] * u_vec_f[i] + C_vec_f[i] * u_vec_f[i+1] - D_vec_f[i];
        }
        check_diff_f[N-1] = A_vec_f[N-1] * u_vec_f[N-2] + B_vec_f[N-1] * u_vec_f[N-1] - D_vec_f[N-1];
    }
    if (to_neumann) {
        check_diff_f[0] = B_vec_f[0] * u_vec_f[0] + (A_vec_f[0] + C_vec_f[0]) * u_vec_f[1] - D_vec_f[0];
        for (i = 1; i < N; i++) {
            check_diff_f[i] = A_vec_f[i] * u_vec_f[i-1] + B_vec_f[i] * u_vec_f[i] + C_vec_f[i] * u_vec_f[i+1] - D_vec_f[i];
        }
        check_diff_f[N] = (A_vec_f[N] + C_vec_f[N]) * u_vec_f[N-1] + B_vec_f[N] * u_vec_f[N] - D_vec_f[N];
    }
}
