#include <stdio.h>
#include <sys/stat.h>
#include <errno.h>
#include <dirent.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

int create_empty_dir(char *parent_directory);

int main()
{
    if (create_empty_dir("./auto_results") != 0) {
        return 1;
    }

    int id = fork();
    if (id == 0) {
        printf("running solver\n");
        execlp("./solver",
               "./solver",
               "./input.txt",
               "./mesh_output.txt",
               "./auto_results",
               "69",
               NULL);

    } else {
        while (wait(NULL) != -1);
        printf("finished\n");
    }

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
