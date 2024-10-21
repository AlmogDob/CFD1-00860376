#include <stdio.h>
#include <sys/stat.h>
#include <errno.h>
#include <dirent.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdarg.h>

int create_empty_dir(char *parent_directory);
void print_command_to_file(FILE *fp, char *program, ...);

int main()
{
    char parent_dir[] = "./auto_results";
    if (create_empty_dir("./auto_results") != 0) {
        return 1;
    }

    char temp_dir[BUFSIZ];
    strncpy(temp_dir, parent_dir, BUFSIZ);
    strncat(temp_dir, "/command_to_run.txt", BUFSIZ/2);
    FILE *fp = fopen(temp_dir, "wt");

    fprintf(fp, "make build_solver\n");

    print_command_to_file(fp,
                          "./solver",
                          "./input.txt",
                          "./mesh_output.txt",
                          "./auto_results",
                          "70",
                          NULL);
    print_command_to_file(fp,
                          "./solver",
                          "./input.txt",
                          "./mesh_output.txt",
                          "./auto_results",
                          "69",
                          NULL);


    // int id = fork();
    // if (id == 0) {
    //     printf("running solver\n");
    //     execlp("./solver",
    //            "./solver",
    //            "./input.txt",
    //            "./mesh_output.txt",
    //            "./auto_results",
    //            "70",
    //            NULL);

    // } else {
    //     while (wait(NULL) != -1);
    //     printf("finished\n");
    // }

    return 0;
}

/* if allready exisest, delet all the files inside 
returns 0 on success
this functin handls the errors so on fail just quit */
int create_empty_dir(char *parent_directory)
{
    char path_to_remove[BUFSIZ];

    if (mkdir(parent_directory, 0777) == -1) {
        if (errno == EEXIST) {
            DIR *dir = opendir(parent_directory);
            if (dir == NULL) {
                fprintf(stderr, "%s:%d: [Error] problem opening '%s': %s\n", __FILE__, __LINE__, parent_directory, strerror(errno));
                return 1;
            }
            struct dirent* entity;
            entity = readdir(dir);
            printf("\n");
            while (entity != NULL) {   /* directory is not empty */
                strncpy(path_to_remove, parent_directory, BUFSIZ);
                strncat(path_to_remove, "/", BUFSIZ/2);
                strncat(path_to_remove, entity->d_name, BUFSIZ);
                printf("%hhd: %s\n", entity->d_type, path_to_remove);
                if (entity->d_type == DT_REG) {
                    if (remove(path_to_remove) != 0) {
                        fprintf(stderr, "%s:%d: [Error] problem removing '%s': %s\n", __FILE__, __LINE__, path_to_remove, strerror(errno));
                        return 1;
                    }
                    printf("remove %s\n", path_to_remove);
                }
                entity = readdir(dir);
            }


            printf("\ndirectory already exist\n\n");

            closedir(dir);

            return 0;
        }

        fprintf(stderr, "%s:%d: [Error] problem making '%s': %s\n", __FILE__, __LINE__, parent_directory, strerror(errno));
        return 1;
    }
    return 0;
}

void print_command_to_file(FILE *fp, char *program, ...)
{
    // printf("%s\n", program);

    // va_list args;
    // va_start(args, program);
    // char *temp_string;
    // while ((temp_string = va_arg(args, char *)) != NULL) {
    //     printf("%s\n", temp_string);
    // }
    // printf("\n");

    // va_end(args);

    fprintf(fp, "%s ", program);
    va_list args;
    va_start(args, program);
    char *temp_string;
    while ((temp_string = va_arg(args, char *)) != NULL) {
        fprintf(fp, "%s ", temp_string);
    }
    va_end(args);
    fprintf(fp, "\n");

}
