
#include <stdio.h>
#include <math.h>
#include "pcgc.h"

void matVecMult(real_t *mat, real_t *vec, real_t *res, int n) 
{
    /* Laço para multiplicar a linha da matriz pelo vetor em forma de coluna */
    for (int i = 0; i < n; ++i) {
        long double sum = 0.0L;
        size_t row_offset = (size_t)i * (size_t)n;
        for (int j = 0; j < n; ++j) {
            sum += (long double)mat[row_offset + j] * (long double)vec[j];
        }
        res[i] = (real_t)sum;
    }
}

real_t dotProduct(real_t *v1, real_t *v2, int n) 
{
    long double sum = 0.0L; 
    /* Laço para carregar o produto escalar em sum */
    for (int i = 0; i < n; ++i) 
        sum += (long double)v1[i] * (long double)v2[i];
    
    return (real_t)sum;
}

void vecAdd(real_t *v1, real_t *v2, real_t scalar, real_t *res, int n) 
{
    /* Laço para carregar a adição vetorial em res */
    for (int i = 0; i < n; ++i) 
        res[i] = v1[i] + scalar * v2[i];
}

real_t maxNormDiff(real_t *v_new, real_t *v_old, int n) 
{
    real_t max_diff = 0.0;
    /* Laço para calcular a diferença das normas e obter a maior */
    for (int i = 0; i < n; ++i) {
        real_t diff = fabs(v_new[i] - v_old[i]);
        if (diff > max_diff) 
            max_diff = diff;
    }
    return max_diff;
}

int solvePCG(real_t *A, real_t *b, real_t *M_inv, real_t *x, int n, long int maxit, 
    real_t epsilon, rtime_t *tempo_iter, real_t *norma_final) 
{
    *tempo_iter = timestamp();

    /* Inicializações dos vetores auxiliares r, z, p, Ap e x_old, e também do iterador k */
    long int k;
    size_t nL = (size_t)n;
    real_t *r = malloc(nL * sizeof(real_t));
    real_t *z = malloc(nL * sizeof(real_t));
    real_t *p = malloc(nL * sizeof(real_t));
    real_t *Ap = malloc(nL * sizeof(real_t));
    real_t *x_old = malloc(nL * sizeof(real_t));
    
    if (!r || !z || !p || !Ap || !x_old) {
        fprintf(stderr, "Erro na alocação de memória na solvePCG\n");
        free(r); free(z); free(p); free(Ap); free(x_old);
        *tempo_iter = timestamp() - *tempo_iter;
        return -1;
    }

    /* Calcula o resíduo inicial r = b - Ax e final r = b + (-1.0 * r) */
    matVecMult(A, x, r, n);
    vecAdd(b, r, -1.0, r, n);

    /* Calcula o pré-condicionador z = M⁻¹r e inicializa p = z */
    matVecMult(M_inv, r, z, n);
    vecAdd(z, z, 0.0, p, n); 

    /* Calcula o produto escalar inicial e zera norma_final para comparação */
    real_t r_dot_z_old = dotProduct(r, z, n);
    *norma_final = 0.0;

    /* Laço principal para iteração do método PCG */
    for (k = 0; k < maxit; ++k) {
        /* Salva x antigo e calcula a norma do erro */
        vecAdd(x, x, 0.0, x_old, n);

        /* Calcula alpha dado por (r'z) / (p'Ap), faz adição vetorial em x = x + alpha * p 
        e subtração em r = r - alpha * Ap */
        matVecMult(A, p, Ap, n);
        real_t alpha = r_dot_z_old / dotProduct(p, Ap, n);
        vecAdd(x, p, alpha, x, n);
        vecAdd(r, Ap, -alpha, r, n);

        /* Calcula a norma do erro para critério de parada */
        real_t actual_norm = maxNormDiff(x, x_old, n);
        *norma_final = actual_norm;

        /* Verifica a norma para critério de parada */
        if (actual_norm < epsilon) {
            free(r); free(z); free(p); free(Ap); free(x_old);
            *tempo_iter = timestamp() - *tempo_iter;
            return k + 1;
        }
       
        /* Calcula z_new = M⁻¹r_new */
        matVecMult(M_inv, r, z, n);
        
        /* Calcula beta = (r_new'z_new) / (r_old'z_old) */
        real_t r_dot_z_new = dotProduct(r, z, n);
        real_t beta = r_dot_z_new / r_dot_z_old;
        r_dot_z_old = r_dot_z_new;

        /* Faz adição vetorial para p_new = z_new + beta * p */
        vecAdd(z, p, beta, p, n);
    }
    
    /* Libera os vetores auxiliares e retorna sem ter convergido */
    free(r); free(z); free(p); free(Ap); // free(x_old);
    *tempo_iter = timestamp() - *tempo_iter;
    fprintf(stderr, "O método não convergiu em %ld iterações.\n", maxit);
    return k + 1;
}