#include "sislin.h"
#include "utils.h"
#include "pcgc.h"

int main () {
    int n, k;
    real_t omega, epsilon;
    long int maxIters;

    /* Lê todos os 5 parâmetros de entrada de uma vez */
    scanf("%d %d %lf %ld %lf", &n, &k, &omega, &maxIters, &epsilon);

    /* Valida os parâmetros lidos */
    if (n <= 10) {
        fprintf(stderr, "n deve ser maior que 10\n");
        return 1;
    }

    if (k <= 1 || k % 2 == 0) {
        fprintf(stderr, "k deve ser maior que 1 e ímpar\n");
        return 1;
    }

    if (omega < -1.0 || omega >= 2.0) {
        fprintf(stderr, "Omega deve estar no intervalo [-1.0, 2.0)\n");
        return 1;
    }

    if (maxIters <= 0) {
        fprintf(stderr, "maxIters deve ser maior que 0\n");
        return 1;
    }

    if (epsilon <= 0.0) {
        fprintf(stderr, "Epsilon deve ser maior que 0\n");
        return 1;
    }

    /* Seed de valores aleatórios */
    srandom(20252);
    printf("\n");

    /* Variáveis para matrizes, vetores, tempos e alocação */
    real_t *A, *b, *ASP, *bsp, *x;
    real_t *D, *L, *U, *M_inv;
    real_t norma_final;
    rtime_t t_pc = 0.0, t_iter = 0.0, t_residuo = 0.0, t_temp;
    size_t nL = (size_t)n;

    /* Gera o sistema k-diagonal e transforma para sistema simétrico positivo */
    criaKDiagonal(n, k, &A, &b);
    bsp = malloc(nL * sizeof(real_t));
    if (!bsp) {
        fprintf(stderr, "Erro de alocação para bsp\n");
        free(A); free(b);
        return 1;
    }

    /* Define a largura de banda para ASP por truncamento */
    int half_a = k / 2;
    int half_asp = 2 * half_a; 
    int k_asp = 2 * half_asp + 1; 

    /* Transforma para sistema simétrico positivo */
    genSimetricaPositiva(A, b, n, k, &ASP, bsp, &t_temp);
    if (!ASP) {
        free(A); free(b); free(bsp);
        return 1;
    }

    /* Aloca o vetor solução */
    x = calloc(nL, sizeof(real_t));
    if (!x) {
        fprintf(stderr, "Erro de alocação para x\n");
        free(A); free(b); free(bsp); free(ASP);
        return 1;
    }

    /* Inicia o tempo para calcular o pré-condicionador (se necessário) e verifica qual 
    pré-condicionador usar a partir do valor de omega */
    t_pc = timestamp();
    if (omega == -1.0) { // Sem pré-condicionador (M = I)
        M_inv = calloc(nL * nL, sizeof(real_t));
        if (!M_inv) {
             fprintf(stderr, "Erro de alocação para M_inv com omega == -1.0\n");
             free(A); free(b); free(bsp); free(ASP); free(x);
             return 1;
        }

        for (int i = 0; i < n; ++i) 
            M_inv[i * n + i] = 1.0;
        
    } else if (omega == 0.0) { // Com pré-condicionador de Jacobi (M = D)
        M_inv = calloc(nL * nL, sizeof(real_t));
        if (!M_inv) {
             fprintf(stderr, "Erro de alocação para M_inv com omega == 0.0\n");
             free(A); free(b); free(bsp); free(ASP); free(x);
             return 1;
        }

        for (int i = 0; i < n; ++i) {
            real_t diag_val = ASP[i * n + i];
            if (fabs(diag_val) < 1e-10) { // Evita divisão por zero
                fprintf(stderr, "Divisão por zero presente na diagonal principal (Jacobi)\n");
                return 1;
            }

            M_inv[i * n + i] = 1.0 / diag_val;
        }

    } else { // Com pré-condicionador Gauss-Seidel ou SSOR (variação do omega dentro de geraPreCond)
        geraDLU(ASP, n, k, &D, &L, &U, &t_temp);
        geraPreCond(D, L, U, omega, n, k, &M_inv, &t_temp);
        if (!M_inv) {
            fprintf(stderr, "Pré-condicionador SSOR/GS não foi gerado\n");
            return 1;
        }

        free(D); free(L); free(U);
    }
    t_pc = timestamp() - t_pc;

    /* Calcula o sistema linear e computa se convergiu ou não */
    int iters = solvePCG(ASP, bsp, M_inv, x, n, maxIters, epsilon, &t_iter, &norma_final);

    /* Caso não tenha convergido */
    if (iters == -1) {
        free(A); free(b);
        free(ASP); free(bsp);
        free(x); free(M_inv);
        return 1;
    }

    /* Calcula o resíduo final e imprime os resultados */
    real_t residuo = calcResiduoSL(ASP, bsp, x, n, k_asp, &t_residuo);
    printf("%d\n", n);
    for (int i = 0; i < n; ++i) 
        printf("%.16g ", x[i]);
    
    printf("\n");
    printf("%.8g\n", norma_final); 
    printf("%.8g\n", residuo);
    printf("%.8g\n", t_pc);
    printf("%.8g\n", t_iter / iters);
    printf("%.8g\n", t_residuo);

    /* Libera a memória restante */
    free(A); free(b);
    free(ASP); free(bsp);
    free(x); free(M_inv);

    return 0;
}