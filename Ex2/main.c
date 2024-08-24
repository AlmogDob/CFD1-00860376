#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <stdlib.h>

#define MAXDIR 1000
#define MAXWORD MAXDIR
#define PI M_PI
#define dprintSTRING(expr) printf(#expr " = %s\n", expr)    /* macro for easy debuging*/
#define dprintINT(expr) printf(#expr " = %d\n", expr)
#define dprintF(expr) printf(#expr " = %g\n", expr)
#define dprintD(expr) printf(#expr " = %g\n", expr)

void read_input(char *input_dir, char *mesh_dir, double *x_vals_mat,
                double *y_vals_mat);
void read_mesh_file(FILE *mesh_fp, double *x_vals_mat,
                   double *y_vals_mat);
void read_mat_from_file(FILE *fp, double *des);
void output_solution(char *dir, double *data);
int offset2d(int i, int j, int ni);
int offset3d(int i, int j, int k, int ni, int nj);
void mat_print(double *data);

/* global variables */
int ni, nj;

int main(int argc, char const *argv[])
{
    /* declerations */
    char input_dir[MAXDIR], mesh_dir[MAXDIR], current_word[MAXWORD];
    double *x_vals_mat, *y_vals_mat;


    /* getting the input directory and mesh directory*/
    if (--argc != 2) {
        fprintf(stderr, "ERROR: not right usage\nUsage: main 'input dir' 'mesh dir'\n");
        return -1;
    }

    strncpy(input_dir, (*(++argv)), MAXDIR);

    if (input_dir[MAXDIR-1] != '\0') {
        fprintf(stderr, "Error: input too long\n");
        return -1;
    }

    strncpy(mesh_dir, (*(++argv)), MAXDIR);

    if (mesh_dir[MAXDIR-1] != '\0') {
        fprintf(stderr, "Error: input too long\n");
        return -1;
    }

    // dprintSTRING(input_dir);
    // dprintSTRING(mesh_dir);

/*------------------------------------------------------------*/

    /* getting ni, nj*/
    FILE *mesh_fp = fopen(mesh_dir, "rt");
    if (!mesh_fp) {
        fprintf(stderr, "Error opening input file: %s\n", strerror(errno));
        exit(1);
    }
    
    while(fscanf(mesh_fp, "%s", current_word) != EOF) {  
        if (!strcmp(current_word, "ni")) {
            fscanf(mesh_fp, "%d", &ni);
        } else if (!strcmp(current_word, "nj")) {
            fscanf(mesh_fp, "%d", &nj);
        }
    }

    x_vals_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (int i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (int j_index = 0; j_index < nj; j_index++) {
            x_vals_mat[offset2d(i_index, j_index, ni)] = 0;
        }
    }
    y_vals_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (int i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (int j_index = 0; j_index < nj; j_index++) {
            y_vals_mat[offset2d(i_index, j_index, ni)] = 0;
        }
    }

/*------------------------------------------------------------*/

    read_input(input_dir, mesh_dir, x_vals_mat, y_vals_mat);

/*------------------------------------------------------------*/

    /* Checking that I got the right input */
    dprintINT(ni);
    dprintINT(nj);
    // printf("x_vals_mat\n");
    // mat_print(x_vals_mat);
    // printf("y_vals_mat\n");
    // mat_print(y_vals_mat);

    printf("--------------------\n");

    /*------------------------------------------------------------*/

    free(x_vals_mat);
    free(y_vals_mat);


    return 0;
}

/* sets 'flags' and variables according to the input file
argument list:
input_dir - the directory of the input file 
mesh_dir - the directory of the mesh file
x_vals_mat - 1D array of the x valuse 
y_vals_mat - 1D array of the y valuse */
void read_input(char *input_dir, char *mesh_dir, double *x_vals_mat,
                double *y_vals_mat)
{
    FILE *input_fp = fopen(input_dir, "rt");
    FILE *mesh_fp = fopen(mesh_dir, "rt");
    char current_word[MAXWORD];

    if (!input_fp) {
        fprintf(stderr, "Error opening input file: %s\n", strerror(errno));
        exit(1);
    }
    if (!mesh_fp) {
        fprintf(stderr, "Error opening input file: %s\n", strerror(errno));
        exit(1);
    }

    read_mesh_file(mesh_fp, x_vals_mat, y_vals_mat);


    fclose(input_fp);
    fclose(mesh_fp);
}

/* reading the mesh file
argument list:
mesh_dir - the directory of the mesh file
x_vals_mat - 1D array of the x valuse 
y_vals_mat - 1D array of the y valuse */
void read_mesh_file(FILE *mesh_fp, double *x_vals_mat,
                   double *y_vals_mat)
{
    char current_word[MAXWORD];

    /* Seting the input varibles according to the mesh file */
    while(fscanf(mesh_fp, "%s", current_word) != EOF) {  
        if (!strcmp(current_word, "x_vals")) {
            read_mat_from_file(mesh_fp, x_vals_mat);
        } else if (!strcmp(current_word, "y_vals")) {
            read_mat_from_file(mesh_fp, y_vals_mat);
        }
    }
}

/* read matrix from file into a 2D matrix (1D array)
fp - file pointer
des - 1D array */
void read_mat_from_file(FILE *fp, double *des)
{
    float temp;

    for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
            fscanf(fp, "%g", &temp);
            des[offset2d(i, j, ni)] = (double)temp;
        }
    }
}

/* converts a 2D index into 1D index
argument list:
i - first direction
j - second direction
ni - first direction size */
int offset2d(int i, int j, int ni)
{
    return j * ni + i;
}

/* converts a 3D index into 1D index
argument list:
i - first direction
j - second direction
k - third direction
ni - first direction size
nj - second direction size */
int offset3d(int i, int j, int k, int ni, int nj)
{
    return (k * nj + j) * ni + i;
}

void mat_print(double *data)
{
    int j_index, i_index;
    for (j_index = 0; j_index < nj; j_index++) {
        for (i_index = 0; i_index < ni; i_index++) {
            printf("%g ", data[offset2d(i_index, j_index, ni)]);
        }
        printf("\n");
    }
    printf("\n");
}