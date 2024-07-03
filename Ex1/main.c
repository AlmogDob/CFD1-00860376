#include <stdio.h>
#define MATRIX_IMPLEMENTATION
#include "Matrix_Double.h"  /* I included this library for debuging */
#include <string.h>
#include <errno.h>
#include <math.h>
#include <assert.h>

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
void interpulat_mat(double *x_vals_mat, char diraction);
double first_deriv(double *mat, char diraction, int i, int j);
double second_deriv(double *mat, char diraction, int i, int j);
void alpha_beta_gama(double *alpha_val_mat, double *beta_val_mat, double *gama_val_mat, double *x_val_mat, double *y_val_mat);

/* Input variables */
double t, delta_x, delta_y, XSF, YSF, x_int = 1;
int i_max, j_max, i_TEL, i_LE, i_TEU, i_min = 0, j_min = 0;

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

    double *alpha_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            alpha_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }

    /*------------------------------------------------------------*/
    
    initialize(x_vals_mat, y_vals_mat);

    double ans = second_deriv(x_vals_mat, 'i', 1, 0);
    dprintD(ans);

    /*------------------------------------------------------------*/

    Mat xmat = {.cols = i_max+1, .rows = j_max+1, .stride = i_max+1, .elements = x_vals_mat}; /* for debuging */
    Mat ymat = {.cols = i_max+1, .rows = j_max+1, .stride = i_max+1, .elements = y_vals_mat}; /* for debuging */
    Mat alphaMat = {.cols = i_max+1, .rows = j_max+1, .stride = i_max+1, .elements = alpha_vals_mat}; 

    // for (i_index = 0; i_index < i_max+1; i_index++) {
    //     alpha_vals_mat[offset2d(i_index, j_max, i_max+1)] = 3;
    // }
    // MAT_PRINT(xmat);

    FILE *fp = fopen("x_mat_output.txt", "wt");
    mat_print_to_file(fp, xmat, "");
    fclose(fp);

    fp = fopen("y_mat_output.txt", "wt");
    mat_print_to_file(fp, ymat, "");
    fclose(fp);

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

/* converts a 2D index into 1D index
argument list: i - x position
j - y position
ni - stride */
int offset2d(int i, int j, int ni)
{
    return j * ni + i;
}

/* set inital valuse of the mash points 
argument list: x_vals_mat - 1D array of the x valuse 
y_vals_mat - 1D array of the y valuse */
void initialize(double *x_vals_mat, double *y_vals_mat)
{
    set_grid_boundaries(x_vals_mat, y_vals_mat);
    interpulat_mat(x_vals_mat, 'j');
    interpulat_mat(y_vals_mat, 'j');
}

/* set the mash boundaries coorditates 
argument list: x_vals_mat - 1D array of the x valuse 
y_vals_mat - 1D array of the y valuse */
void set_grid_boundaries(double *x_vals_mat, double *y_vals_mat)
{
    int i_index, i_min = 0, j_index, j_min = 0, index = 0, num_points_befor_circle,
    num_of_outer_segments, num_of_top_outer_segments;

    double x, x_i_minos_2, x_i_minos_1, y_j_minos_2, y_j_minos_1, y_imax_jmax, x_imax_jmax,
    delta_theta, R, theta = 0, length_top_outer, segment_length, current_x_val = 0;

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
    /* Outer boundary */
    y_imax_jmax = y_vals_mat[offset2d(i_max, j_max, i_max+1)];
    x_imax_jmax = x_vals_mat[offset2d(i_max, j_max, i_max+1)];
    R = y_imax_jmax;

    num_of_outer_segments = i_max;
    num_of_top_outer_segments = num_of_outer_segments/2;
    length_top_outer = x_imax_jmax + 0.5*PI*R; /* length of stright part and quarter of the circle */
    segment_length = length_top_outer/num_of_top_outer_segments;

    /* the stright line part */
    for (num_points_befor_circle = 0;
         num_points_befor_circle < num_of_top_outer_segments + 1;
         num_points_befor_circle++) {
            current_x_val = x_imax_jmax - num_points_befor_circle*segment_length; 
            if (current_x_val < 0) {
                break;
            }
            x_vals_mat[offset2d(i_max-num_points_befor_circle, j_max, i_max+1)] = current_x_val;
            y_vals_mat[offset2d(i_max-num_points_befor_circle, j_max, i_max+1)] = y_imax_jmax;
    }

    theta = PI/2 + atan(x_vals_mat[offset2d(i_max-num_points_befor_circle+1, j_max, i_max+1)] / R);
    delta_theta = theta / (num_of_top_outer_segments - num_points_befor_circle + 1);

    /* the quarter circle part */
    for (index = 0; index < num_of_top_outer_segments - num_points_befor_circle + 1; index++) {
        x_vals_mat[offset2d(i_max/2 + index, j_max, i_max+1)] = -R*cos(delta_theta*index);
        y_vals_mat[offset2d(i_max/2 + index, j_max, i_max+1)] = R*sin(delta_theta*index);
    }

    /* coping to the bottom side */
    for (index = 1; index < i_max/2 + 1; index++) {
        x_vals_mat[offset2d(i_max/2 - index, j_max, i_max+1)] = x_vals_mat[offset2d(i_max/2 + index, j_max, i_max+1)];
        y_vals_mat[offset2d(i_max/2 - index, j_max, i_max+1)] = -y_vals_mat[offset2d(i_max/2 + index, j_max, i_max+1)];
    }
}

/* returns the shape of the airfoil as a function of x
argument list: x - x positon 
side - the uppers side are the low side of the airfoil */
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

/* fill the matrix in a way of interpulation the boundarys
argument list: mat - 1D array of valuse
diraction - i direction or j direction */
void interpulat_mat(double *mat, char diraction)
{
    int i, j; 
    double max, min;

    assert(diraction == 'j' || diraction == 'i'); /* different directions are not implemented */

    if (diraction == 'j') {
        for (i = 1; i < i_max; i++) {
            for (j = 1; j < j_max; j++) {
                max = mat[offset2d(i, j_max, i_max+1)];
                min = mat[offset2d(i, j_min, i_max+1)];
                mat[offset2d(i, j, i_max+1)] = (max - min)/(j_max) * j + min;   /* liniar interpulation */
            }
        }
    }
    if (diraction == 'i') {
        for (j = 1; j < j_max; j++) {
            for (i = 1; i < i_max; i++) {
                max = mat[offset2d(i_max, j, i_max+1)];
                min = mat[offset2d(i_min, j, i_max+1)];
                mat[offset2d(i, j, i_max+1)] = (max - min)/(i_max) * i + min;   /* liniar interpulation */
            }
        }
    }
}

double first_deriv(double *mat, char diraction, int i, int j)
{
    if (diraction == 'j') {
        assert(j != j_min && j != j_max);
        return (mat[offset2d(i, j+1, i_max+1)] - mat[offset2d(i, j-1, i_max+1)]) / (2); /* second order first derivitive */
    }
    if (diraction == 'i') {
        assert(i != i_min && i != i_max);
        return (mat[offset2d(i+1, j, i_max+1)] - mat[offset2d(i-1, j, i_max+1)]) / (2); /* second order first derivitive */
    }
}

double second_deriv(double *mat, char diraction, int i, int j)
{
    if (diraction == 'j') {
        assert(j != j_min && j != j_max);
        return (mat[offset2d(i, j+1, i_max+1)] -2*mat[offset2d(i, j, i_max)] + mat[offset2d(i, j-1, i_max+1)]) / (1); /* second order second derivitive */
    }
    if (diraction == 'i') {
        assert(i != i_min && i != i_max);
        return (mat[offset2d(i+1, j, i_max+1)] -2*mat[offset2d(i, j, i_max)] + mat[offset2d(i-1, j, i_max+1)]) / (1); /* second order second derivitive */
    }
}

