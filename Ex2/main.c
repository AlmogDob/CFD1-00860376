#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <stdlib.h>

#define MAXDIR 1000
#define MAXWORD MAXDIR
#define PI M_PI
#define dprintSTRING(expr) do{printf("%d: ", __LINE__); printf(#expr " = %s\n", expr);} while(0)   /* macro for easy debuging*/
#define dprintINT(expr) do{printf("%d: ", __LINE__); printf(#expr " = %d\n", expr);} while(0)     /* macro for easy debuging*/
#define dprintF(expr) do{printf("%d: ", __LINE__); printf(#expr " = %g\n", expr);} while(0)     /* macro for easy debuging*/
#define dprintD(expr) do{printf("%d: ", __LINE__); printf(#expr " = %g\n", expr);} while(0)     /* macro for easy debuging*/

void read_input(char *input_dir, char *mesh_dir, double *x_vals_mat,
                double *y_vals_mat);
void read_mesh_file(FILE *mesh_fp, double *x_vals_mat,
                   double *y_vals_mat);
void read_input_file(FILE *fp);
void read_mat_from_file(FILE *fp, double *des);
void output_solution(char *dir, double *data);
int offset2d(int i, int j, int ni);
int offset3d(int i, int j, int k, int ni, int nj);
void print_mat2D(double *data);
void print_layer_of_mat3D(double *data, int layer);
double first_deriv(double *mat, char diraction, int i, int j);
double second_deriv(double *mat, char diraction, int i, int j);
double calculate_one_over_jacobian_at_a_point(double *x_vals_mat,
                                              double *y_vals_mat, int i,
                                              int j);
void contravariant_velocities(double *U, double *V, double *x_vals_mat,
                              double *y_vals_mat, double *Q, int i, int j);
void calculate_u_and_v(double *u, double *v, double *Q, int i, int j);
double calculate_energy(void);
void calculate_E_hat_at_a_point(double *E0, double *E1, double *E2,
                                double *E3, double *x_vals_mat,
                                double *y_vals_mat, double *Q, int i,
                                int j);
void calculate_F_hat_at_a_point(double *F0, double *F1, double *F2,
                                double *E3, double *x_vals_mat,
                                double *y_vals_mat, double *Q, int i,
                                int j);                               
void initialize_flow_field(double *Q);
void RHS(double *S, double *W, double *Q, double *x_vals_mat, double *y_vals_mat);
void advance_Q(double *next_Q, double *current_Q ,double *S, double *x_vals_mat,
               double *y_vals_mat);

/* global variables */
int ni, nj, max_ni_nj, i_TEL, i_LE, i_TEU, j_TEL, j_LE, j_TEU;
double Mach, angle_of_attack_deg, angle_of_attack_rad, density,
environment_pressure, delta_t, Gamma;

int main(int argc, char const *argv[])
{
/* declerations */
    char input_dir[MAXDIR], mesh_dir[MAXDIR], current_word[MAXWORD];
    double *x_vals_mat, *y_vals_mat, *J_vals_mat, *current_Q, *next_Q, *S, *W;

/* getting the input directory and mesh directory*/
    if (--argc != 2) {
        fprintf(stderr, "%s:%d: [Error] not right usage... Usage: main 'input dir' 'mesh dir'\n", __FILE__, __LINE__);
        return -1;
    }

    strncpy(input_dir, (*(++argv)), MAXDIR);

    if (input_dir[MAXDIR-1] != '\0') {
        fprintf(stderr, "%s:%d: [Error] input too long\n", __FILE__, __LINE__);
        return -1;
    }

    strncpy(mesh_dir, (*(++argv)), MAXDIR);

    if (mesh_dir[MAXDIR-1] != '\0') {
        fprintf(stderr, "%s:%d: [Error] input too long\n", __FILE__, __LINE__);
        return -1;
    }

    // dprintSTRING(input_dir);
    // dprintSTRING(mesh_dir);

/*------------------------------------------------------------*/

/* getting ni, nj*/
    FILE *mesh_fp = fopen(mesh_dir, "rt");
    if (!mesh_fp) {
        fprintf(stderr, "%s:%d: [Error] opening input file: %s\n",__FILE__, __LINE__, strerror(errno));
        exit(1);
    }
    
    while(fscanf(mesh_fp, "%s", current_word) != EOF) {  
        if (!strcmp(current_word, "ni")) {
            fscanf(mesh_fp, "%d", &ni);
        } else if (!strcmp(current_word, "nj")) {
            fscanf(mesh_fp, "%d", &nj);
        }
    }
    fclose(mesh_fp);
    max_ni_nj = (int)fmax((double)ni, (double)nj);

/*------------------------------------------------------------*/
/* allocating the matrices */
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
    J_vals_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (int i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (int j_index = 0; j_index < nj; j_index++) {
            J_vals_mat[offset2d(i_index, j_index, ni)] = 0;
        }
    }
    W = (double *)malloc(sizeof(double) * max_ni_nj * 4);
    for (int i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        for (int j_index = 0; j_index < 4; j_index++) {
            W[offset2d(i_index, j_index, ni)] = 0;
        }
    }
    current_Q = (double *)malloc(sizeof(double) * ni * nj * 4);
    for (int i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (int j_index = 0; j_index < nj; j_index++) {
            for (int k_index = 0; k_index < 4; k_index++) {
                current_Q[offset3d(i_index, j_index, k_index, ni, nj)] = 0;
            }
        }
    }
    next_Q = (double *)malloc(sizeof(double) * ni * nj * 4);
    for (int i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (int j_index = 0; j_index < nj; j_index++) {
            for (int k_index = 0; k_index < 4; k_index++) {
                next_Q[offset3d(i_index, j_index, k_index, ni, nj)] = 0;
            }
        }
    }
    S = (double *)malloc(sizeof(double) * ni * nj * 4);
    for (int i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (int j_index = 0; j_index < nj; j_index++) {
            for (int k_index = 0; k_index < 4; k_index++) {
                S[offset3d(i_index, j_index, k_index, ni, nj)] = 0;
            }
        }
    }

/*------------------------------------------------------------*/

    read_input(input_dir, mesh_dir, x_vals_mat, y_vals_mat);

/* Checking the input */
    printf("--------------------\n");
    dprintINT(ni);
    dprintINT(nj);
    // printf("x_vals_mat\n");
    // print_mat2D(x_vals_mat);
    // printf("y_vals_mat\n");
    // print_mat2D(y_vals_mat);
    dprintINT(i_TEL);
    dprintINT(i_LE);
    dprintINT(i_TEU);
    dprintINT(j_TEL);
    dprintINT(j_LE);
    dprintINT(j_TEU);
    dprintD(Mach);
    dprintD(angle_of_attack_deg);
    dprintD(density);
    dprintD(environment_pressure);
    dprintD(delta_t);
    dprintD(Gamma);
    dprintINT(max_ni_nj);

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            J_vals_mat[offset2d(i, j, ni)] = 1.0 / calculate_one_over_jacobian_at_a_point(x_vals_mat, y_vals_mat, i, j);
            // J_vals_mat[offset2d(i, j, ni)] = i + j;
        }
    }
    print_mat2D(J_vals_mat);

    printf("--------------------\n");

/*------------------------------------------------------------*/

    initialize_flow_field(current_Q);

    
    RHS(S, W, current_Q, x_vals_mat, y_vals_mat);
    advance_Q(next_Q, current_Q, S, x_vals_mat, y_vals_mat);
    
    // int layer = 3;
    // print_layer_of_mat3D(current_Q, layer);
    // print_layer_of_mat3D(next_Q, layer);

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

    if (!input_fp) {
        fprintf(stderr, "%s:%d: [Error] failed opening input file: %s\n", __FILE__, __LINE__, strerror(errno));
        exit(1);
    }
    if (!mesh_fp) {
        fprintf(stderr, "%s:%d: [Error] failed opening mesh file: %s\n", __FILE__, __LINE__, strerror(errno));
        exit(1);
    }

    read_mesh_file(mesh_fp, x_vals_mat, y_vals_mat);
    read_input_file(input_fp);

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

/* read input parameters from input file
fp - file pointer */
void read_input_file(FILE *fp)
{
    char current_word[MAXWORD];
    float temp;

    while(fscanf(fp, "%s", current_word) != EOF) {  
        if (!strcmp(current_word, "Mach")) {
            fscanf(fp, "%g", &temp);
            Mach = (double)temp;
        } else if (!strcmp(current_word, "angle_of_attack_deg")) {
            fscanf(fp, "%g", &temp);
            angle_of_attack_deg = (double)temp;
            angle_of_attack_rad = PI*angle_of_attack_deg/180.0;
        } else if (!strcmp(current_word, "density")) {
            fscanf(fp, "%g", &temp);
            density = (double)temp;
        } else if (!strcmp(current_word, "environment_pressure")) {
            fscanf(fp, "%g", &temp);
            environment_pressure = (double)temp;
        } else if (!strcmp(current_word, "delta_t")) {
            fscanf(fp, "%g", &temp);
            delta_t = (double)temp;
        } else if (!strcmp(current_word, "Gamma")) {
            fscanf(fp, "%g", &temp);
            Gamma = (double)temp;
        } else if (!strcmp(current_word, "i_TEL")) {
            fscanf(fp, "%d", &i_TEL);
        } else if (!strcmp(current_word, "i_LE")) {
            fscanf(fp, "%d", &i_LE);
        } else if (!strcmp(current_word, "i_TEU")) {
            fscanf(fp, "%d", &i_TEU);
        } else if (!strcmp(current_word, "j_TEL")) {
            fscanf(fp, "%d", &j_TEL);
        } else if (!strcmp(current_word, "j_LE")) {
            fscanf(fp, "%d", &j_LE);
        } else if (!strcmp(current_word, "j_TEU")) {
            fscanf(fp, "%d", &j_TEU);
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

void print_mat2D(double *data)
{
    int j_index, i_index;

    for (j_index = nj - 1; j_index >= 0; j_index--) {
        for (i_index = 0; i_index < ni; i_index++) {
            printf("%g ", data[offset2d(i_index, j_index, ni)]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_layer_of_mat3D(double *data, int layer)
{
    int j_index, i_index;

    for (j_index = nj - 1; j_index >= 0; j_index--) {
        for (i_index = 0; i_index < ni; i_index++) {
            printf("%g ", data[offset3d(i_index, j_index, layer, ni, nj)]);
        }
        printf("\n");
    }
    printf("\n");
}

/* return the second order first derivitive from the valuse in the mat matrix
argument list:
mat - 1D array of valuse
diraction - i or j
i, j - the points coordinates */
double first_deriv(double *mat, char diraction, int i, int j)
{
    int j_min = 0, j_max = nj-1, i_min = 0, i_max = ni-1;

    if (diraction == 'j') {
        if (j == j_min || j == j_max) {
            return 0;
        }
        return (mat[offset2d(i, j+1, i_max+1)] - mat[offset2d(i, j-1, i_max+1)]) / (2); /* second order first derivitive */
    }
    if (diraction == 'i') {
        if (i == i_min || i == i_max) {
            return 0;
        }
        return (mat[offset2d(i+1, j, i_max+1)] - mat[offset2d(i-1, j, i_max+1)]) / (2); /* second order first derivitive */
    }
    return NAN;
}

/* return the second order second derivitive from the valuse in the mat matrix
argument list:
mat - 1D array of valuse
diraction - i or j
i, j - the points coordinates */
double second_deriv(double *mat, char diraction, int i, int j)
{
    int j_min = 0, j_max = nj-1, i_min = 0, i_max = ni-1;

    if (diraction == 'j') {
        if (j == j_min || j == j_max) {
            return 0;
        }
        return (mat[offset2d(i, j+1, i_max+1)] -2*mat[offset2d(i, j, i_max+1)] + mat[offset2d(i, j-1, i_max+1)]) / (1); /* second order second derivitive */
    }
    if (diraction == 'i') {
        if (i == i_min || i == i_max) {
            return 0;
        }
        return (mat[offset2d(i+1, j, i_max+1)] -2*mat[offset2d(i, j, i_max+1)] + mat[offset2d(i-1, j, i_max+1)]) / (1); /* second order second derivitive */
    }
    return NAN;
}

/* calculating the jacobian in a single point
argument list:
x_vals_mat - 1D array of the x valuse 
y_vals_mat - 1D array of the y valuse 
i, j - the points coordinates */
double calculate_one_over_jacobian_at_a_point(double *x_vals_mat,
                                              double *y_vals_mat, int i,
                                              int j)
{
    double dx_dxi, dx_deta, dy_dxi, dy_deta;

    dx_dxi = first_deriv(x_vals_mat, 'i', i, j);
    dx_deta = first_deriv(x_vals_mat, 'j', i, j);
    dy_dxi = first_deriv(y_vals_mat, 'i', i, j);
    dy_deta = first_deriv(y_vals_mat, 'j', i, j);

    return (dx_dxi*dy_deta - dy_dxi*dx_deta);
}

/* calculating the contravariant velocities in a single point
argument list:
U - the contravariant velocity in the xi direction
V - the contravariant velocity in the eta direction
x_vals_mat - 1D array of the x valuse 
y_vals_mat - 1D array of the y valuse 
Q - 1D matrix of the flow field
i, j - the points coordinates */
void contravariant_velocities(double *U, double *V, double *x_vals_mat,
                              double *y_vals_mat, double *Q, int i, int j)
{
    double J, dx_dxi, dx_deta, dy_dxi, dy_deta, dxi_dx, dxi_dy, deta_dx,
    deta_dy, u, v;

    J = 1.0 / calculate_one_over_jacobian_at_a_point(x_vals_mat, y_vals_mat, i, j);
    dx_dxi = first_deriv(x_vals_mat, 'i', i, j);
    dx_deta = first_deriv(x_vals_mat, 'j', i, j);
    dy_dxi = first_deriv(y_vals_mat, 'i', i, j);
    dy_deta = first_deriv(y_vals_mat, 'j', i, j);

    dxi_dx  =   J * dy_deta;
    dxi_dy  = - J * dx_deta;
    deta_dx = - J * dy_dxi;
    deta_dy =   J * dx_dxi;

    calculate_u_and_v(&u, &v, Q, i, j);

    *U = dxi_dx  * u + dxi_dy  * v;
    *V = deta_dx * u + deta_dy * v;
}

void calculate_u_and_v(double *u, double *v, double *Q, int i, int j)
{
    *u = Q[offset3d(i, j, 1, ni, nj)] / Q[offset3d(i, j, 0, ni, nj)]; /* rho*u / rho */
    *v = Q[offset3d(i, j, 2, ni, nj)] / Q[offset3d(i, j, 0, ni, nj)]; /* rho*v / rho */
}

double calculate_energy(void)
{
    double p, speed_of_sound, velocity, internal_energy, energy;

    p = environment_pressure;
    speed_of_sound = sqrt(Gamma * p / density);

    velocity = Mach * speed_of_sound;
    internal_energy = p / (Gamma - 1) / density;
    energy = density * internal_energy + density * velocity * velocity / 2;
    return energy;
}

void calculate_E_hat_at_a_point(double *E0, double *E1, double *E2,
                                double *E3, double *x_vals_mat,
                                double *y_vals_mat, double *Q, int i,
                                int j)
{
    double u, v, U, V, one_over_J, dx_deta, dy_deta, dxi_dx, dxi_dy,
    energy, p;

    calculate_u_and_v(&u, &v, Q, i, j);
    contravariant_velocities(&U, &V, x_vals_mat, y_vals_mat, Q, i, j);
    one_over_J = calculate_one_over_jacobian_at_a_point(x_vals_mat, y_vals_mat, i, j);
    dx_deta = first_deriv(x_vals_mat, 'j', i, j);
    dy_deta = first_deriv(y_vals_mat, 'j', i, j);
    dxi_dx  =   1 / one_over_J * dy_deta;
    dxi_dy  = - 1 / one_over_J * dx_deta;
    energy = calculate_energy();
    p = environment_pressure;

    if (!one_over_J) {
        *E0 = 0;
        *E1 = 0;
        *E2 = 0;
        *E3 = 0;
    } else {
        *E0 = one_over_J * density * U;
        *E1 = one_over_J * (density * u * U + dxi_dx * p); 
        *E2 = one_over_J * (density * v * U + dxi_dy * p);
        *E3 = one_over_J * (energy + p) * U;
    }

}

void calculate_F_hat_at_a_point(double *F0, double *F1, double *F2,
                                double *F3, double *x_vals_mat,
                                double *y_vals_mat, double *Q, int i,
                                int j)
{
    double u, v, U, V, one_over_J, dx_dxi, dy_dxi, deta_dx, deta_dy,
    energy, p;

    calculate_u_and_v(&u, &v, Q, i, j);
    contravariant_velocities(&U, &V, x_vals_mat, y_vals_mat, Q, i, j);  
    one_over_J = calculate_one_over_jacobian_at_a_point(x_vals_mat, y_vals_mat, i, j);
    dx_dxi = first_deriv(x_vals_mat, 'i', i, j);
    dy_dxi = first_deriv(y_vals_mat, 'i', i, j);
    deta_dx = - 1 / one_over_J * dy_dxi;
    deta_dy =   1 / one_over_J * dx_dxi;
    energy = calculate_energy();
    p = environment_pressure;

    if (!one_over_J) {
        *F0 = 0;
        *F1 = 0;
        *F2 = 0;
        *F3 = 0;
    } else {
        *F0 = one_over_J * V;
        *F1 = one_over_J * (density * u * V + deta_dx * p);
        *F2 = one_over_J * (density * v * V + deta_dy * p);
        *F3 = one_over_J * (energy + p) * V;
    }

}

void initialize_flow_field(double *Q)
{
    double u, v, energy, p, speed_of_sound, velocity;

    p = environment_pressure;
    speed_of_sound = sqrt(Gamma * p / density);

    velocity = Mach * speed_of_sound;
    u = velocity * cos(angle_of_attack_rad);
    v = velocity * sin(angle_of_attack_rad);

    energy = calculate_energy();

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            Q[offset3d(i, j, 0, ni, nj)] = density;
            Q[offset3d(i, j, 1, ni, nj)] = density * u;
            Q[offset3d(i, j, 2, ni, nj)] = density * v;
            Q[offset3d(i, j, 3, ni, nj)] = energy;
        }
    }
}

void RHS(double *S, double *W, double *Q, double *x_vals_mat, double *y_vals_mat)
{
    int i, j, k;

    /* zeroing S and W*/
    for (i = 0; i < ni; i++) {   
        for (j = 0; j < nj; j++) {
            for (k = 0; k < 4; k++) {
                S[offset3d(i, j, k, ni, nj)] = 0;
            }
        }
    }
    for (i = 0; i < max_ni_nj; i++) {
        for (j = 0; j < 4; j++) {
            W[offset2d(i, j, max_ni_nj)] = 0;
        }
    }

    /* xi direction (constant j) */
    for (j = 1; j < nj - 1; j++) {
        for (i = 0; i < ni; i++) {
            calculate_E_hat_at_a_point(&W[offset2d(i, 0, max_ni_nj)],
                                       &W[offset2d(i, 1, max_ni_nj)],
                                       &W[offset2d(i, 2, max_ni_nj)],
                                       &W[offset2d(i, 3, max_ni_nj)],
                                       x_vals_mat, y_vals_mat, Q, i, j);
        }

        for (i = 1; i < ni - 1; i++) {
            for (k = 0; k < 4; k++) {
                S[offset3d(i, j, k, ni, nj)] += -delta_t * 0.5 * (W[offset2d(i+1, k, max_ni_nj)] - W[offset2d(i-1, k, max_ni_nj)]);
            }
        }
    } 

    /* eta direction (constant i) */
    for (i = 1; i < ni - 1; i++) {
        for (j = 0; j < nj; j++) {
            calculate_F_hat_at_a_point(&W[offset2d(j, 0, max_ni_nj)],
                                       &W[offset2d(j, 1, max_ni_nj)],
                                       &W[offset2d(j, 2, max_ni_nj)],
                                       &W[offset2d(j, 3, max_ni_nj)],
                                       x_vals_mat, y_vals_mat, Q, i, j);
        }

        for (j = 1; j < nj - 1; j++) {
            for (k = 0; k < 4; k++) {
                S[offset3d(i, j, k, ni, nj)] += -delta_t * 0.5 * (W[offset2d(j+1, k, max_ni_nj)] - W[offset2d(j-1, k, max_ni_nj)]);
            }
        }
    }

    // for (j = 0; j < 4; j++) {
    //     for (i = 0; i < max_ni_nj; i++) {
    //         printf("%g ", W[offset2d(i, 3-j, max_ni_nj)]);
    //     }
    //     printf("\n");
    // }

    // print_layer_of_mat3D(S, 2);

}

void advance_Q(double *next_Q, double *current_Q ,double *S, double *x_vals_mat,
               double *y_vals_mat)
{
    double J;

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            J = 1.0 / calculate_one_over_jacobian_at_a_point(x_vals_mat, y_vals_mat, i, j);
            for (int k = 0; k < 4; k++ ) {
                int index = offset3d(i, j, k, ni, nj);
                if (S[index]) {
                    next_Q[index] = current_Q[index] + S[index] * J;
                } else {
                    next_Q[index] = current_Q[index];
                }
            }
        }
    }
}
