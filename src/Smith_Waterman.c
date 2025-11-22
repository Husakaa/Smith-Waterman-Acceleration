#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>


#define MAX_SEQ_LEN 20000
#define GAP -2
#define DIAGONAL 1
#define IZQUIERDA 2
#define ARRIBA 3
#define CERO 0

int caracter(char base) {
    switch(base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

int max4(int a, int b, int c, int d) {
    int max = a;
    if (b > max) max = b;
    if (c > max) max = c;
    if (d > max) max = d;
    return max;
}

void leer_fasta(const char *ruta, char *secuencia) {
    FILE *f = fopen(ruta, "r");
    if (!f) {
        fprintf(stderr, "Error al abrir el archivo: %s\n", ruta);
        exit(EXIT_FAILURE);
    }

    char linea[1024];
    secuencia[0] = '\0';

    while (fgets(linea, sizeof(linea), f)) {
        if (linea[0] == '>') continue;
        linea[strcspn(linea, "\r\n")] = 0;
        if (strlen(secuencia) + strlen(linea) < MAX_SEQ_LEN)
            strcat(secuencia, linea);
        else {
            fprintf(stderr, "Secuencia demasiado larga para buffer\n");
            exit(EXIT_FAILURE);
        }
    }

    fclose(f);
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Uso: %s <fichero_fasta1> <fichero_fasta2>\n", argv[0]);
        return EXIT_FAILURE;
    }


    char secuencia1[MAX_SEQ_LEN];
    char secuencia2[MAX_SEQ_LEN];

    leer_fasta(argv[1], secuencia1);
    leer_fasta(argv[2], secuencia2);

    int m = strlen(secuencia1);
    int n = strlen(secuencia2);

    printf("La longitud de %s es %d\n", argv[1], m);
    printf("La longitud de %s es %d\n", argv[2], n);

    int similitud[4][4] = {
        {3, -1, 1, -1}, // A
        {-1, 3, -1, 1}, // C
        {1, -1, 3, -1}, // G
        {-1, 1, -1, 3}  // T
    };


    // Reservar matrices dinámicas
    int *H = malloc((m + 1) * (n + 1) * sizeof(int));
    int *traceback = malloc((m + 1) * (n + 1) * sizeof(int));
    if (!H || !traceback) {
        fprintf(stderr, "No se pudo asignar memoria para las matrices\n");
        return EXIT_FAILURE;
    }

    // Inicialización
    for (int i = 0; i <= m; i++) {
        for (int j = 0; j <= n; j++) {
            H[i * (n + 1) + j] = 0;
            traceback[i * (n + 1) + j] = 0;
        }
    }

    double t_real_start = omp_get_wtime();  // Tiempo real total
    clock_t t_cpu_start = clock();          // Tiempo de CPU total
    // Rellenar matrices
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int s = 0;
            int ind1 = caracter(secuencia1[i - 1]);
            int ind2 = caracter(secuencia2[j - 1]);
            if (ind1 >= 0 && ind2 >= 0)
                s = similitud[ind1][ind2];

            int diag = H[(i - 1) * (n + 1) + (j - 1)] + s;
            int up = H[(i - 1) * (n + 1) + j] + GAP;
            int left = H[i * (n + 1) + (j - 1)] + GAP;
            int val = max4(0, diag, up, left);
            H[i * (n + 1) + j] = val;

            if (val == 0)
                traceback[i * (n + 1) + j] = CERO;
            else if (val == diag)
                traceback[i * (n + 1) + j] = DIAGONAL;
            else if (val == up)
                traceback[i * (n + 1) + j] = ARRIBA;
            else
                traceback[i * (n + 1) + j] = IZQUIERDA;
        }
    }

    clock_t t_cpu_end = clock();
    double t_real_end = omp_get_wtime();

    double cpu_seconds = (double)(t_cpu_end - t_cpu_start) / CLOCKS_PER_SEC;
    double real_seconds = t_real_end - t_real_start;

    // Buscar el máximo
    int max_i = 0, max_j = 0, max_score = 0;
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int val = H[i * (n + 1) + j];
            if (val > max_score) {
                max_score = val;
                max_i = i;
                max_j = j;
            }
        }
    }

    // Reconstrucción
    char alineacion1[MAX_SEQ_LEN];
    char alineacion2[MAX_SEQ_LEN];
    int i = max_i, j = max_j, pos = 0;

    while (traceback[i * (n + 1) + j] != CERO) {
        if (traceback[i * (n + 1) + j] == DIAGONAL) {
            alineacion1[pos] = secuencia1[i - 1];
            alineacion2[pos] = secuencia2[j - 1];
            i--; j--;
        } else if (traceback[i * (n + 1) + j] == ARRIBA) {
            alineacion1[pos] = secuencia1[i - 1];
            alineacion2[pos] = '-';
            i--;
        } else {
            alineacion1[pos] = '-';
            alineacion2[pos] = secuencia2[j - 1];
            j--;
        }
        pos++;
    }

    alineacion1[pos] = '\0';
    alineacion2[pos] = '\0';

    // Invertir
    for (int k = 0; k < pos / 2; k++) {
        char tmp = alineacion1[k];
        alineacion1[k] = alineacion1[pos - 1 - k];
        alineacion1[pos - 1 - k] = tmp;

        tmp = alineacion2[k];
        alineacion2[k] = alineacion2[pos - 1 - k];
        alineacion2[pos - 1 - k] = tmp;
    }

   

    printf("Tiempo real total: %.3f segundos\n", real_seconds);
    printf("Tiempo de CPU total: %.3f segundos\n", cpu_seconds);
    printf("Score: %d", max_score);

    // Liberar memoria
    free(H);
    free(traceback);

    return 0;
}
