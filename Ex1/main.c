#include <stdio.h>
#define MATRIX_IMPLEMENTATION
#include "Matrix_Double.c"
#include <string.h>

#define MAXDIR 100
#define dprintSTRING(expr) printf(#expr " = %s\n", expr)    /* macro for easy debuging*/
#define dprintINT(expr) printf(#expr " = %d\n", expr)
#define dprintF(expr) printf(#expr " = %g\n", expr)
#define dprintD(expr) printf(#expr " = %g\n", expr)

int main(int argc, char const *argv[])
{
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

    dprintSTRING(input_dir);
    dprintSTRING(output_dir);


    return 0;
}
