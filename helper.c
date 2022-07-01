#include "helper.h"
#include <stdio.h>

void show_matrix(unsigned char *matrix, int rows, int cols, char *message)
{
    printf("-- %s ----------\n", message);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%d, ", matrix[i * cols + j]);
        }
        printf("\n");
    }
}

void show_frags(unsigned char **frags, int rows, char *message)
{
    printf("-- %s ----------\n", message);
    for (int i = 0; i < rows; i++)
    {
        printf("%d: ", i);
        printf("(%s)\n", frags[i]);
    }
}