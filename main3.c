#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/*
El formato del archivo matrix.txt debe ser el siguiente:

n
a11 a12 a13 ... a1n
a21 a22 a23 ... a2n
..
an1 an2 an3 ... ann
*/

void print_matrix(double *A, int  tamMatrix) {
    for (int i = 0; i < tamMatrix; i++) {
        for (int j = 0; j < tamMatrix; j++) {
            printf("%f ", A[i * tamMatrix + j]);
        }
        printf("\n");
    }
}

int compare(const void *a, const void *b) {
    return (*(double *)a > *(double *)b) - (*(double *)a < *(double *)b);
}

int main(int argc, char **argv) {
    int rank, nProc, tamMatrix;
    double *A = NULL;
    double *sub_A;
    double *sums = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);

    if (rank == 0) {
        // Leer la matriz A desde el archivo
        FILE *file = fopen("matrix.txt", "r");
        if (file == NULL) {
            fprintf(stderr, "Error al abrir el archivo matrix.txt\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        fscanf(file, "%d", &tamMatrix); // Leer el tamaño de la matriz
        A = (double *)malloc(tamMatrix * tamMatrix * sizeof(double));
        
        for (int i = 0; i < tamMatrix; i++) {
            for (int j = 0; j < tamMatrix; j++) {
                fscanf(file, "%lf", &A[i * tamMatrix + j]); // Leer elementos
            }
        }
        fclose(file);

        // printf("Matriz A inicial:\n");
        // print_matrix(A, tamMatrix);

        sums = (double *)malloc(tamMatrix * sizeof(double)); // Array para guardar sumas
        
        // Calcular la distribución de filas entre procesos
        int rows_per_process = tamMatrix / nProc;
        int extra_rows = tamMatrix % nProc;
        int *sendcounts = (int *)malloc(nProc * sizeof(int));
        int *displs = (int *)malloc(nProc * sizeof(int));

        // Asignar sendcounts y displs para cada proceso
        int offset = 0;
        for (int i = 0; i < nProc; i++) {
            if (i < extra_rows) {
                sendcounts[i] = (rows_per_process + 1) * tamMatrix;
            } else {
                sendcounts[i] = rows_per_process * tamMatrix;
            }
            displs[i] = offset;
            offset += sendcounts[i];
        }

        /*
        printf("--Antes de scatterv--");
        printf("rank: %d\n", rank);
        printf("n: %d\n", tamMatrix);
        printf("nProc: %d\n", nProc);
        printf("rows_per_process: %d\n", rows_per_process);
        printf("extra_rows: %d\n", extra_rows);
        printf("sendcounts: ");
        for (int i = 0; i < nProc; i++) {
            printf("%d ", sendcounts[i]);
        }
        printf("\n");
        printf("displs: ");
        for (int i = 0; i < nProc; i++) {
            printf("%d ", displs[i]);
        }
        printf("\n");
        */
    }
    
    // Distribuir filas de A a cada proceso usando MPI_Scatterv
    MPI_Scatterv(A, sendcounts, displs, MPI_DOUBLE, sub_A, local_rows * tamMatrix, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Calcular la suma de cada fila local
    double *local_sums = (double *)malloc(local_rows * sizeof(double));
    for (int i = 0; i < local_rows; i++) {
        local_sums[i] = 0.0;
        for (int j = 0; j < n; j++) {
            local_sums[i] += sub_A[i * tamMatrix  + j];
        }
    }

    // Recolectar todas las sumas en el proceso raíz (rank 0)
    int *recvcounts = NULL;
    int *displs_sums = NULL;

    if (rank == 0) {
        recvcounts = (int *)malloc(nProc * sizeof(int));
        displs_sums = (int *)malloc(nProc * sizeof(int));

        offset = 0;
        for (int i = 0; i < nProc; i++) {
            if (i < extra_rows) {
                recvcounts[i] = rows_per_process + 1; // Si hay filas extra, asigna una fila más
            } else {
                recvcounts[i] = rows_per_process; // De lo contrario, asigna el número estándar de filas
            }
            displs_sums[i] = offset; // Asigna el desplazamiento actual
            offset += recvcounts[i]; // Actualiza el desplazamiento
        }
    }

    // Ahora, usando recvcounts y displs_sums solo si rank es 0
    MPI_Gatherv(local_sums, local_rows, MPI_DOUBLE, sums, recvcounts, displs_sums, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    if (rank == 0) {
        // Ordenar las filas según las sumas
        double *sorted_A = (double *)malloc(tamMatrix  * tamMatrix * sizeof(double));
        int *order = (int *)malloc(tamMatrix * sizeof(int));

        // Asignar las sumas a un arreglo de estructuras para ordenación
        for (int i = 0; i < tamMatrix; i++) {
            order[i] = i;
        }

        // Ordenar las filas según las sumas
        qsort(order, tamMatrix, sizeof(int), compare);

        // Reorganizar la matriz A según el orden
        for (int i = 0; i < tamMatrix; i++) {
            int idx = order[i];
            for (int j = 0; j < tamMatrix; j++) {
                sorted_A[i * tamMatrix + j] = A[idx * tamMatrix + j];
            }
        }

        // Escribir la matriz ordenada en un archivo
        FILE *output_file = fopen("sorted_matrix.txt", "w");
        for (int i = 0; i < tamMatrix; i++) {
            for (int j = 0; j < tamMatrix; j++) {
                fprintf(output_file, "%f ", sorted_A[i * tamMatrix + j]);
            }
            fprintf(output_file, "\n");
        }
        fclose(output_file);

        // Limpiar memoria
        free(A);
        free(sums);
        free(sorted_A);
        free(order);
        free(recvcounts);
        free(displs_sums);
    }

    // Limpiar memoria
    free(sub_A);
    free(local_sums);
    free(sendcounts);
    free(displs);
    */
    MPI_Finalize();
    return 0;
}
