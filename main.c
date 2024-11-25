#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>    
#include <sys/time.h> 
/*
INPUT del programa:
    - Ruta del archivo donde se almacena la matriz a ordenar.
    - Numero de procesos a utilizar.


El formato del archivo de entrada debe ser el siguiente:

n
a11 a12 a13 ... a1n
a21 a22 a23 ... a2n
..
an1 an2 an3 ... ann

La matriz ordenada se escribe en el archivo sorted_matrix.txt en el directorio
actual. En el archivo sorted_matrix.txt, todos los elementos son doubles y se 
representan con todos sus decimales.
*/

/*
    Función para imprimir una matriz
    INPUT: 
        - A: matriz a imprimir
        - tamMatrix: tamaño de la matriz
*/
void print_matrix(double *A, int tamMatrix) {
    for (int i = 0; i < tamMatrix; i++) {
        for (int j = 0; j < tamMatrix; j++) {
            printf("%f ", A[i * tamMatrix + j]);
        }
        printf("\n");
    }
}

static double *sums;

/*
    Función de comparación para 'qsort'
    INPUT: 
        - a: primer elemento a comparar
        - b: segundo elemento a comparar
*/
int compare_sums(const void *a, const void *b) {
    int index_a = *(const int *)a;
    int index_b = *(const int *)b;
    return (sums[index_a] > sums[index_b]) - (sums[index_a] < sums[index_b]);
}


int main(int argc, char **argv) {
    int rank, nProc, tamMatrix;
    if (argc != 3) {
        fprintf(stderr, "Numero de argumentos insuficientes\n");
        fprintf(stderr, "\tArgumento 1: Numero de procesos a utilizar\n");
        fprintf(stderr, "\tArgumento 2: Ruta del fichero con la matriz a ordenar\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    nProc = atoi(argv[1]);
    char* filename = argv[2];

    double *A = NULL;
    double *sub_A = NULL;
    double *A_cyclic = NULL;
    
    
    // Variables necesarias para la distribución
    int *sendcounts = NULL; // sendcounts: cantidad de datos que se envían a cada proceso
    int *displs = NULL; // displs: desplazamiento de los datos enviados a cada proceso
    int extra_rows = 0; // Filas adicionales que se distribuyen entre los procesos

    // Variables para medir tiempos de ejecucion
	struct timeval start;
	struct timeval start2;
	struct timeval end;

    // Inicializar MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    
    if (rank == 0){
        // Leer la matriz desde el archivo matrix.txt
        FILE *file = fopen(filename, "r");
        if (file == NULL) {
            fprintf(stderr, "Error al abrir el archivo %s\n", filename);
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

        // Array para guardar sumas
        sums = (double *)malloc((tamMatrix) * sizeof(double));  // +1 porque el proceso 0 también envia datos (0)
        // Inicializamos sums a -1
        for (int i = 0; i < tamMatrix; i++) {
            sums[i] = -1;
        }

        // Calcular la distribución de filas entre procesos
        int rows_per_process;
        rows_per_process = tamMatrix / (nProc - 1); // Quitamos al proceso 0
        extra_rows = tamMatrix % (nProc - 1);

        
        sendcounts = (int *)malloc(nProc * sizeof(int));
        displs = (int *)malloc(nProc * sizeof(int));

        // Crear una matriz temporal para reorganizar las filas de manera cíclica
        A_cyclic = (double *)malloc(tamMatrix * tamMatrix * sizeof(double));
        // Reordenar la matriz A de manera cíclica
        int row_index = 0;
        for (int proc = 1; proc < nProc; proc++) {
            for (int i = 0; i < tamMatrix; i++) {
                if (i % (nProc - 1) == (proc - 1)) { // Verifica si la fila i le toca al proceso actual
                    // Copiar la fila i de A a la posición actual en A_cyclic
                    for (int j = 0; j < tamMatrix; j++) {
                        A_cyclic[row_index * tamMatrix + j] = A[i * tamMatrix + j];
                    }
                    row_index++; // Avanzar al siguiente índice de fila en A_cyclic
                }
            }
        }

        

        #ifdef DEBUG
        printf("Matriz reordenada;");
        for (int i = 0; i < tamMatrix*tamMatrix; i++){
            printf("%lf ", A_cyclic[i]);
            printf("\n");
        }
        printf("\n");
        #endif

        // Marca el inicio de la medicion del tiempo
        gettimeofday(&start, NULL);
        gettimeofday(&start2, NULL);

        // Asignar sendcounts y displs para cada proceso
        int offset = 0;
        for (int i = 0; i < nProc; i++) {
            if (i == 0) {
                sendcounts[i] = 0; // El proceso 0 no recibe datos
                displs[0] = 0; // Asignar el desplazamiento actual
            } 
            else{
                if (i < extra_rows + 1) {
                sendcounts[i] = (rows_per_process + 1) * tamMatrix; // Un proceso recibe una fila más
                } else {
                    sendcounts[i] = rows_per_process * tamMatrix; // Procesos reciben filas estándar
                }
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


    double *local_sums;
    sub_A = (double *)malloc(sendcounts[rank] * sizeof(double));    
    MPI_Scatterv(A_cyclic, sendcounts, displs, MPI_DOUBLE, sub_A, sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    #ifdef DEBUG
    if (rank == 0) {
        
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
        
    }
    printf("--- Proceso %d, sendcount %d - displs %d ---\n", rank, *sendcounts, displs[rank]);
    
    // Imprimimos la matriz que le llega a cada proceso
    printf("\tProceso %d, sub_A:", rank);
    for (int i = 0; i < sendcounts[rank]; i++) {
        printf(" %f", sub_A[i]);
    }
    printf("\n");
    #endif

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

    #ifdef DEBUG
    printf("\tProceso %d, suma local:", rank);
    for (int i = 0; i < local_rows; i++) {
        printf(" %f", local_sums[i]);
    }
    printf("\n");
    #endif
    
    ///////// LAS SUMAS FUNCIONAN
    
    // Recolectar todas las sumas en el proceso raíz (rank 0) utilizando MPI_Gatherç
    int rows_per_process[nProc];
    rows_per_process[0] = 0;
    for (int i = 1; i < nProc; i++) {
        rows_per_process[i] = sendcounts[i] / tamMatrix;
    }

    int indiceSums[nProc]; // El 0 no se usa (envia 0 filas)
    // int seSumo = 0; // Variable para saber si se sumo una fila extra en el anterior indiceSums
    if (rank == 0){
        indiceSums[0] = 0;
        for (int i = 1; i < nProc; i++) {
            indiceSums[i] = indiceSums[i-1] + rows_per_process[i-1];
        }
    }
    

    MPI_Gatherv(local_sums, rows_per_process[rank], MPI_DOUBLE, sums, rows_per_process, indiceSums, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #ifdef DEBUG
    if (rank == 0){
        printf("\t\t\t Parametros de Gatherv\n");
        printf("\t\t\t - rows_per_process: ");
        for (int i = 0; i < nProc; i++) {
            printf("%d ", rows_per_process[i]);
        }
        printf("\n");
        printf("\t\t\t - indiceSums: ");
        for (int i = 0; i < nProc + 1; i++) {
            printf("%d ", indiceSums[i]);
        }
        printf("\n");

        printf("\t\t Proceso 0 - Sumas totales: ");
        for (int i = 0; i < tamMatrix; i++) {
            printf("%f ", sums[i]);
        }
        printf("\n");
    }
    #endif

    // SUMAS TOTALES FUNCIONAN
    
    
    if (rank == 0) {
        // Ordenar las filas según las sumas
        double *sorted_A = (double *)malloc(tamMatrix * tamMatrix * sizeof(double));
        int *order = (int *)malloc(tamMatrix * sizeof(int));

        // Inicializar el arreglo 'order' para el índice de las filas
        for (int i = 0; i < tamMatrix; i++) {
            order[i] = i;
        }

        

        // Ordenar las filas según las sumas usando 'qsort'
        qsort(order, tamMatrix, sizeof(int), compare_sums);

        #ifdef DEBUG
        printf("order\n");
        for (int i = 0; i < tamMatrix; i++) {
            printf("%d ", order[i]);
        }
        printf("\n");
        #endif

        // Reorganizar la matriz A según el orden de 'sums'
        for (int i = 0; i < tamMatrix; i++) {
            int idx = order[i];
            for (int j = 0; j < tamMatrix; j++) {
                sorted_A[i * tamMatrix + j] = A_cyclic[idx * tamMatrix + j];
            }
        }


        // Marca el fin de la medicion del tiempo
        gettimeofday(&end, NULL);
        // Calcula el tiempo de sobrecarga antes del inicio real del calculo
        double overhead = (start2.tv_sec - start.tv_sec) + (start2.tv_usec - start.tv_usec) / 1e6;
        // Calcula el tiempo total del calculo excluyendo la sobrecarga
        double time = (end.tv_sec - start2.tv_sec) + (end.tv_usec - start2.tv_usec) / 1e6 - overhead;

        #ifdef DEBUG
        // Escribir la matriz ordenada en un archivo
        FILE *output_file = fopen("sorted_matrix.txt", "w");
        for (int i = 0; i < tamMatrix; i++) {
            for (int j = 0; j < tamMatrix; j++) {
                fprintf(output_file, "%f ", sorted_A[i * tamMatrix + j]);
            }
            fprintf(output_file, "\n");
        }
        fclose(output_file);
        #endif
        

        // Escribimos resultados de tiempos:
        FILE *tiempos = fopen("data.csv", "a");
        fprintf(tiempos, "%d;%d;%lf\n", tamMatrix, nProc, time);
        fclose(tiempos);


        free(sorted_A);
        free(order);
        free(A);
        free(sums);
        free(A_cyclic);
        sorted_A = NULL;
        order = NULL;
        A = NULL;
        sums = NULL;
        A_cyclic = NULL;
    }

    free(sub_A);
    free(local_sums);
    free(sendcounts);
    free(displs);
    sub_A = NULL;
    local_sums = NULL;
    sendcounts = NULL;
    displs = NULL;

    
    
    
    MPI_Finalize();
    return 0;
}
