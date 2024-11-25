#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

// Función para intercambiar dos filas de la matriz
void swap_rows(int **matrix, int row1, int row2, int size) {
    for (int i = 0; i < size; i++) {
        int temp = matrix[row1][i];
        matrix[row1][i] = matrix[row2][i];
        matrix[row2][i] = temp;
    }
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Uso: %s <tamaño de la matriz>\n", argv[0]);
        return 1;
    }

    int size = atoi(argv[1]);
    if (size <= 0) {
        fprintf(stderr, "El tamaño de la matriz debe ser un número positivo.\n");
        return 1;
    }

    // Reservar memoria para la matriz
    int **matrix = (int **)malloc(size * sizeof(int *));
    for (int i = 0; i < size; i++) {
        matrix[i] = (int *)malloc(size * sizeof(int));
    }

    // Rellenar la matriz con números doubles incrementales
    int value = 1;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            matrix[i][j] = value++;
        }
    }

    // Inicializar el generador de números aleatorios
    srand(time(NULL));

    // Reordenar las filas de la matriz en orden aleatorio
    for (int i = 0; i < size; i++) {
        int random_row = rand() % size;
        swap_rows(matrix, i, random_row, size);
    }

    char nombreArchivo[255];
    strcpy(nombreArchivo, "./matrices/"); // Copiar la ruta base
    strcat(nombreArchivo, "matrix");     // Concatenar el prefijo del archivo
    strcat(nombreArchivo, argv[1]);      // Concatenar el argumento proporcionado
    strcat(nombreArchivo, ".txt");       // Concatenar la extensión

    // Abrir el archivo para escritura
    FILE *file = fopen(nombreArchivo, "w");

    if (file == NULL) {
        fprintf(stderr, "Error al abrir el archivo.\n");
        for (int i = 0; i < size; i++) {
            free(matrix[i]);
        }
        free(matrix);
        return 1;
    }

    fprintf(file, "%d\n",size);
    int limite = size - 1; 
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (j != limite){
                fprintf(file, "%d ", matrix[i][j]); // Imprimir los números con precisión decimal
            }
            else{
                fprintf(file, "%d", matrix[i][j]); // Imprimir los números con precisión decimal
            }
        }
        fprintf(file, "\n");
    }

    fclose(file);

    // Liberar memoria
    for (int i = 0; i < size; i++) {
        free(matrix[i]);
    }
    free(matrix);


    printf("\t -- Matriz de %dx%d creada --\n", size, size);
    return 0;
}