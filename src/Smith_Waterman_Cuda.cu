#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>

#define MAX_SEQ_LEN 20000
#define GAP -2

// Tabla de similitud (A, C, G, T) en memoria constante
__constant__ int d_similitud[16];

// Función para convertir base a índice
__host__ __device__ int caracter(char base) {
    switch (base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

// Leer archivo FASTA
void leer_fasta(const char *ruta, char *secuencia) {
    FILE *f = fopen(ruta, "r");
    if (!f) {
        fprintf(stderr, "Error al abrir %s\n", ruta);
        exit(EXIT_FAILURE);
    }
    char linea[1024];
    secuencia[0] = '\0';
    while (fgets(linea, sizeof(linea), f)) {
        if (linea[0] == '>') continue;
        linea[strcspn(linea, "\r\n")] = '\0';
        strcat(secuencia, linea);
    }
    fclose(f);
}

// Kernel CUDA: calcular diagonal diag de la matriz H
__global__ void smith_waterman_kernel(
    char *seq1, char *seq2, int *H, int m, int n, int diag, int start_i, int end_i) 
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i = start_i + idx;
    if (i <= end_i) {
        int j = diag - i + 1;
        int ind1 = caracter(seq1[i - 1]);
        int ind2 = caracter(seq2[j - 1]);
        int s = (ind1 >= 0 && ind2 >= 0) ? d_similitud[ind1 * 4 + ind2] : 0;

        int idx_cur = i * (n + 1) + j;
        int idx_diag = (i - 1) * (n + 1) + (j - 1);
        int idx_up = (i - 1) * (n + 1) + j;
        int idx_left = i * (n + 1) + (j - 1);

        int score_diag = H[idx_diag] + s;
        int score_up = H[idx_up] + GAP;
        int score_left = H[idx_left] + GAP;

        int max_score = max(0, max(score_diag, max(score_up, score_left)));
        H[idx_cur] = max_score;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        printf("Uso: %s <fasta1> <fasta2> <blockSize>\n", argv[0]);
        return 1;
    }

    // Medición tiempo total (CPU)
    clock_t start_total = clock();

    char *secuencia1 = (char *)malloc(MAX_SEQ_LEN);
    char *secuencia2 = (char *)malloc(MAX_SEQ_LEN);

    leer_fasta(argv[1], secuencia1);
    leer_fasta(argv[2], secuencia2);
    int threadsPerBlock = atoi(argv[3]);
    if (threadsPerBlock % 32 != 0) {
        printf("Block size must be multiple of 32");
        return 1;
    }

    int m = strlen(secuencia1);
    int n = strlen(secuencia2);

    printf("Longitud %s: %d\n", argv[1], m);
    printf("Longitud %s: %d\n", argv[2], n);

    size_t size = (m + 1) * (n + 1) * sizeof(int);
    int *H_host = (int *)calloc((m + 1) * (n + 1), sizeof(int));

    char *d_seq1, *d_seq2;
    int *d_H;

    // Medición tiempo GPU + transferencias con cudaEvent_t
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    cudaMalloc(&d_seq1, m);
    cudaMalloc(&d_seq2, n);
    cudaMalloc(&d_H, size);

    
    cudaMemcpy(d_seq1, secuencia1, m, cudaMemcpyHostToDevice);
    cudaMemcpy(d_seq2, secuencia2, n, cudaMemcpyHostToDevice);
    cudaMemset(d_H, 0, size);

    int h_similitud[16] = {
         3, -1,  1, -1,
        -1,  3, -1,  1,
         1, -1,  3, -1,
        -1,  1, -1,  3
    };
    cudaMemcpyToSymbol(d_similitud, h_similitud, sizeof(int) * 16);


    // Recorremos todas las diagonales (desde 1 hasta m+n-1)
    for (int diag = 1; diag <= m + n - 1; diag++) {
        int start_i = diag - n + 1;
        if (start_i < 1) start_i = 1;
        int end_i = diag;
        if (end_i > m) end_i = m;

        int num_elements = end_i - start_i + 1;
        if (num_elements <= 0) continue;

        int blocksPerGrid = (num_elements + threadsPerBlock - 1) / threadsPerBlock;

        smith_waterman_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_seq1, d_seq2, d_H, m, n, diag, start_i, end_i);
        cudaDeviceSynchronize();
    }

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);

    float ms_gpu;
    cudaEventElapsedTime(&ms_gpu, start, stop);

    cudaMemcpy(H_host, d_H, size, cudaMemcpyDeviceToHost);

    // Buscar el score máximo
    int max_score = 0;
    for (int i = 0; i <= m; i++) {
        for (int j = 0; j <= n; j++) {
            int val = H_host[i * (n + 1) + j];
            if (val > max_score) max_score = val;
        }
    }

    // Fin medición tiempo total (CPU)
    clock_t end_total = clock();
    double tiempo_total = (double)(end_total - start_total) / CLOCKS_PER_SEC;
    printf("Tiempo total programa: %.3f segundos\n", tiempo_total);
    printf("Tiempo GPU + transferencias: %.3f ms\n", ms_gpu);
    printf("Puntuacion maxima: %d\n", max_score);

    // Liberar recursos
    cudaFree(d_seq1);
    cudaFree(d_seq2);
    cudaFree(d_H);
    free(secuencia1);
    free(secuencia2);
    free(H_host);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    return 0;
}
