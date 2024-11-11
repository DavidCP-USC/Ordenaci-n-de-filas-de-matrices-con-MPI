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

void print_matrix(double *A, int tamMatrix) {
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
    nProc = atoi(argv[1]);
    double *A = NULL;
    double *sub_A = NULL; // Inicializar sub_A
    double *sums = NULL;
    
    // Variables necesarias para la distribución
    int *sendcounts = NULL; // sendcounts: cantidad de datos que se envían a cada proceso
    int *displs = NULL; // displs: desplazamiento de los datos enviados a cada proceso

    // Inicializar MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    
    if (rank == 0){
        // Leer la matriz desde el archivo matrix.txt
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

        sums = (double *)malloc(tamMatrix * sizeof(double)); // Array para guardar sumasmalloc(tamMatrix * sizeof(double)); // Array para guardar sumas

        // Calcular la distribución de filas entre procesos
        int rows_per_process, extra_rows;
        rows_per_process = tamMatrix / nProc;
        extra_rows = tamMatrix % nProc;
    
        sendcounts = (int *)malloc(nProc * sizeof(int));
        displs = (int *)malloc(nProc * sizeof(int));
    

        // Asignar sendcounts y displs para cada proceso
        int offset = 0;
        for (int i = 0; i < nProc; i++) {
            if (i < extra_rows) {
                sendcounts[i] = (rows_per_process + 1) * tamMatrix; // Un proceso recibe una fila más
            } else {
                sendcounts[i] = rows_per_process * tamMatrix; // Procesos reciben filas estándar
            }
            displs[i] = offset; // Asignar el desplazamiento actual
            offset += sendcounts[i]; // Actualizar el desplazamiento
        }
    }
    else{
        sendcounts = (int *)malloc(nProc * sizeof(int));
        displs = (int *)malloc(nProc * sizeof(int));
    }

    // Transmitir sendcounts y displs a todos los procesos
    MPI_Bcast(&tamMatrix, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // SE LE ESTA PASANDO A TODOS LA INFO DE TODOS, CAMBIAR
    // TODO: Cambiar a MPI_Scatter
    MPI_Bcast(sendcounts, nProc, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(displs, nProc, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Mostrar detalles para depuración (opcional)
    printf("--Antes de scatterv--\n");
    printf("\trank: %d\n", rank);
    printf("\tn: %d\n", tamMatrix);
    printf("\tnProc: %d\n", nProc);
    printf("\tsendcounts: ");
    for (int i = 0; i < nProc; i++) {
        printf("%d ", sendcounts[i]);
    }
    printf("\n");
    printf("\tdispls: ");
    for (int i = 0; i < nProc; i++) {
        printf("%d ", displs[i]);
    }
    printf("\n");

    double *local_sums;
    sub_A = (double *)malloc(sendcounts[rank] * sizeof(double));    
    MPI_Scatterv(A, sendcounts, displs, MPI_DOUBLE, sub_A, sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        /*
        // Distribuir filas de A a cada proceso usando MPI_Scatterv
        printf("----- Proceso %d -----\n", rank);
        printf("Paramtros del scatterv\n");
        printf("\tA:\n");
        print_matrix(A, tamMatrix);
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
    else{
        printf("--- Proceso %d, sendcount %d - displs %d ---\n", rank, *sendcounts, displs[rank]);
        
        // Imprimimos la matriz que le llega a cada proceso
        printf("\tProceso %d, sub_A:", rank);
        for (int i = 0; i < sendcounts[rank]; i++) {
            printf(" %f", sub_A[i]);
        }
        printf("\n");

        // Definimos local_rows basado en el tamaño que recibimos
        int local_rows = sendcounts[rank] / tamMatrix; // Calculamos cuántas filas locales tiene este proceso

        // Calcular la suma de cada fila local
        local_sums = (double *)malloc(local_rows * sizeof(double));
        for (int i = 0; i < local_rows; i++) {
            local_sums[i] = 0.0;
            for (int j = 0; j < tamMatrix; j++) {
                local_sums[i] += sub_A[i * tamMatrix + j]; // Sumar los elementos de la fila
            }
        }
        printf("\tProceso %d, suma local:", rank);
        for (int i = 0; i < local_rows; i++) {
            printf(" %f", local_sums[i]);
        }
        printf("\n");
    }
    ///////// LAS SUMAS FUNCIONAN
/*
    // Recolectar todas las sumas en el proceso raíz (rank 0) utilizando MPI_Gather
    MPI_Gather(local_sums, rows_per_process, MPI_DOUBLE, sums, rows_per_process, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0){
        printf("Sumas totales: ");
        for (int i = 0; i < tamMatrix; i++) {
            printf("%f ", sums[i]);
        }
        printf("\n");
    }



    if (rank == 0) {
        // Recolectar todas las sumas en el proceso raíz (rank 0)
        int *recvcounts = NULL;
        int *displs_sums = NULL;

        recvcounts = (int *)malloc(nProc * sizeof(int));
        displs_sums = (int *)malloc(nProc * sizeof(int));

        int offset = 0;
        for (int i = 0; i < nProc; i++) {
            if (i < extra_rows) {
                recvcounts[i] = (rows_per_process + 1); // Si hay filas extra, asigna una fila más
            } else {
                recvcounts[i] = rows_per_process; // De lo contrario, asigna el número estándar de filas
            }
            displs_sums[i] = offset; // Asigna el desplazamiento actual
            offset += recvcounts[i]; // Actualiza el desplazamiento
        }
    }
    /*
    if (rank == 0) {
        // Ordenar las filas según las sumas
        double *sorted_A = (double *)malloc(tamMatrix * tamMatrix * sizeof(double));
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
    if (rank == 0) {
        free(sendcounts);
        free(displs);
    } else {
        free(sendcounts);
        free(displs);
    }
    */
    
    MPI_Finalize();
    return 0;
}
