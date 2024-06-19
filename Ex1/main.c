#include <stdio.h>
#define MATRIX_IMPLEMENTATION
#include "Matrix_Double.h"  /* I included this library for debuging */
#include <string.h>
#include <errno.h>
#include <math.h>

#define MAXDIR 100
#define MAXWORD 100
#define PI 3.14159265359
#define dprintSTRING(expr) printf(#expr " = %s\n", expr)    /* macro for easy debuging*/
#define dprintINT(expr) printf(#expr " = %d\n", expr)
#define dprintF(expr) printf(#expr " = %g\n", expr)
#define dprintD(expr) printf(#expr " = %g\n", expr)

void read_input(char *dir);
void output_solution_d(char *dir, double *data);
int offset2d(int i, int j, int ni);
void initialize(double *x_vals_mat, double *y_vals_mat);
void set_grid_boundaries(double *x_vals_mat, double *y_vals_mat);
double airfoil(double x, char side);

/* Input variables */
double t, delta_x, delta_y, XSF, YSF, x_int = 1;
int i_max, j_max, i_TEL, i_LE, i_TEU;

int main(int argc, char const *argv[])
{
    /* decleraitons */
    char input_dir[MAXDIR], output_dir[MAXDIR];
    int i_index, j_index;

    /*------------------------------------------------------------*/

    /* Geting the input and output directories */
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

    /*------------------------------------------------------------*/

    /* Checking that I got the right input */
    dprintD(t);
    dprintINT(i_max);
    dprintINT(j_max);
    dprintINT(i_TEL);
    dprintINT(i_LE);
    dprintINT(i_TEU);
    dprintD(delta_x);
    dprintD(delta_y);
    dprintD(XSF);
    dprintD(YSF);
    printf("--------------------\n");

    /*------------------------------------------------------------*/
    
    /* Memory allocation */
    double *x_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            x_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }

    double *y_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            x_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }

    Mat xmat = {.cols = i_max+1, .rows = j_max+1, .stride = i_max+1, .elements = x_vals_mat}; /* for debuging */
    Mat ymat = {.cols = i_max+1, .rows = j_max+1, .stride = i_max+1, .elements = y_vals_mat}; /* for debuging */


    /*------------------------------------------------------------*/
    
    initialize(x_vals_mat, y_vals_mat);
    /*------------------------------------------------------------*/
    MAT_PRINT(xmat);
    MAT_PRINT(ymat);

    return 0;
}

/* sets 'flags' and variables according to the input file
argument list: dir - the directory of the input file */
void read_input(char *dir)
{
    FILE *fp = fopen(dir, "rt");
    char current_word[MAXWORD];
    float temp;

    if (!fp) {
        fprintf(stderr, "Error opening file: %s\n", strerror(errno));
        exit(1);
    }

    /* Seting the input varibles according to the input file */
    while(fscanf(fp, "%s", current_word) != EOF) {  
        if (!strcmp(current_word, "t")) {
            fscanf(fp, "%g", &temp); /* fscanf returns a float, so I cast it to double */
            t = (double)temp;
        } else if (!strcmp(current_word, "i_max")) {
            fscanf(fp, "%d", &i_max);
        } else if (!strcmp(current_word, "j_max")) {
            fscanf(fp, "%d ", &j_max);
        } else if (!strcmp(current_word, "i_TEL")) {
            fscanf(fp, "%d", &i_TEL);
        } else if (!strcmp(current_word, "i_LE")) {
            fscanf(fp, "%d", &i_LE);
        } else if (!strcmp(current_word, "i_TEU")) {
            fscanf(fp, "%d", &i_TEU);
        } else if (!strcmp(current_word, "delta_y")) {
            fscanf(fp, "%g", &temp); /* fscanf returns a float, so I cast it to double */
            delta_y = (double)temp;
        } else if (!strcmp(current_word, "XSF")) {
            fscanf(fp, "%g", &temp); /* fscanf returns a float, so I cast it to double */
            XSF = (double)temp;
        } else if (!strcmp(current_word, "YSF")) {
            fscanf(fp, "%g", &temp); /* fscanf returns a float, so I cast it to double */
            YSF = (double)temp;
        } else if (!strcmp(current_word, "x_int")) {
            fscanf(fp, "%g", &temp); /* fscanf returns a float, so I cast it to double */
            x_int = (double)temp;
        }
    }

    delta_x = 1.0/(i_LE - i_TEL);

    fclose(fp);
}

/* output data; double version
argument list: dir - the directory of the output file.
data - the solution vector */
void output_solution_d(char *dir, double *data)
{
    FILE *fp = fopen(dir, "wt");
    
    data = (void *)data;

    fclose(fp);
}

int offset2d(int i, int j, int ni)
{
    return j * ni + i;
}

void initialize(double *x_vals_mat, double *y_vals_mat)
{
    set_grid_boundaries(x_vals_mat, y_vals_mat);
}

void set_grid_boundaries(double *x_vals_mat, double *y_vals_mat)
{
    int i_index, i_min = 0, j_index, j_min = 0;
    double x, x_i_minos_2, x_i_minos_1, y_j_minos_2, y_j_minos_1;

    /* setting the boundary according to the exercie */
    /* Eq 6 */
    for (i_index = i_TEL, j_index = j_min; i_index < i_TEU+1; i_index++) { 
        x = 1 - cos(0.5*PI*(i_LE-i_index)*delta_x);
        x_vals_mat[offset2d(i_index, j_index, i_max+1)] = x;
        if (i_index < i_LE) {
            y_vals_mat[offset2d(i_index, j_index, i_max+1)] = airfoil(x, 'l');
        } else if (i_index >= i_LE) {
            y_vals_mat[offset2d(i_index, j_index, i_max+1)] = airfoil(x, 'u');
        }
    }
    /* Eq 7 */
    for (i_index = i_TEU + 1, j_index = j_min; i_index < i_max+1; i_index++) {
        x_i_minos_1 = x_vals_mat[offset2d(i_index-1, j_index, i_max+1)];  
        x_i_minos_2 = x_vals_mat[offset2d(i_index-2, j_index, i_max+1)];  
        x_vals_mat[offset2d(i_index, j_index, i_max+1)] = x_i_minos_1 + (x_i_minos_1 - x_i_minos_2) * XSF;
    }
    for (i_index = i_min, j_index = j_min; i_index < i_TEL; i_index++) {
        x_vals_mat[offset2d(i_index, j_index, i_max+1)] = x_vals_mat[offset2d(i_max-i_index, j_index, i_max+1)];
    }
    /* Eq 8 */
    y_vals_mat[offset2d(i_max, j_min+1, i_max+1)] = delta_y;
    for (i_index = i_max, j_index = j_min+2; j_index < j_max+1; j_index++) {
        y_j_minos_1 = y_vals_mat[offset2d(i_index, j_index-1, i_max+1)];  
        y_j_minos_2 = y_vals_mat[offset2d(i_index, j_index-2, i_max+1)];  
        y_vals_mat[offset2d(i_max, j_index, i_max+1)] = y_j_minos_1 + (y_j_minos_1 - y_j_minos_2) * YSF;
    }
    for (i_index = i_max, j_index = j_min+1; j_index < j_max+1; j_index++) {
        x_vals_mat[offset2d(i_max, j_index, i_max+1)] = x_vals_mat[offset2d(i_max, i_min, i_max+1)];
        y_vals_mat[offset2d(i_min, j_index, i_max+1)] = -y_vals_mat[offset2d(i_max, j_index, i_max+1)];
        x_vals_mat[offset2d(i_min, j_index, i_max+1)] = x_vals_mat[offset2d(i_min, j_min, i_max+1)];
    }
}

double airfoil(double x, char side)
{
    if (side == 'u') {
        return 5 * t * (0.2969*sqrt(x_int*x) - 0.1260*(x_int*x)
               - 0.3516*pow(x_int*x,2) + 0.2843*pow(x_int*x,3)
               - 0.1015*pow(x_int*x,4));
    } else if (side == 'l') {
        return - 5 * t * (0.2969*sqrt(x_int*x) - 0.1260*(x_int*x)
               - 0.3516*pow(x_int*x,2) + 0.2843*pow(x_int*x,3)
               - 0.1015*pow(x_int*x,4));
    }
    fprintf(stderr, "Error with the airfoil usage\n");
    exit(1);
}
