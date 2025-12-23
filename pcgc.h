#ifndef __PCGC_H__
#define __PCGC_H__

#include "utils.h"

/** 
 * Multiplica uma matriz quadrada n x n por um vetor n, armazenando o resultado em outro vetor n
 * @param mat ponteiro para a matriz de coeficientes (armazenada em formato linear).
 * @param vec ponteiro para o vetor a ser multiplicado.
 * @param res ponteiro para o vetor resultado da multiplicação.
 * @param n dimensão da matriz e dos vetores.
 */
void matVecMult(real_t *mat, real_t *vec, real_t *res, int n);

/**
 * Calcula o produto escalar entre dois vetores de dimensão n
 * @param v1 ponteiro para o primeiro vetor.
 * @param v2 ponteiro para o segundo vetor.
 * @param n dimensão dos vetores.
 * @return real_t resultado do produto escalar.
 */
real_t dotProduct(real_t *v1, real_t *v2, int n);

/**
 * Realiza a operação de adição vetorial com escalares: res = v1 + scalar * v2
 * @param v1 ponteiro para o primeiro vetor.
 * @param v2 ponteiro para o segundo vetor.
 * @param scalar escalar multiplicador do segundo vetor.
 * @param res ponteiro para o vetor resultado.
 * @param n dimensão dos vetores.
 */
void vecAdd(real_t *v1, real_t *v2, real_t scalar, real_t *res, int n);

/**
 * Calcula a norma máxima da diferença entre dois vetores de dimensão n
 * @param v_new ponteiro para o vetor novo.
 * @param v_old ponteiro para o vetor antigo.
 * @param n dimensão dos vetores.
 * @return real_t norma máxima da diferença.
 */
real_t maxNormDiff(real_t *v_new, real_t *v_old, int n);

/**
 * Resolve um sistema linear Ax=b usando o método dos Gradientes Conjugados Pré-condicionados
 * @param A matriz de coeficientes (simétrica e positiva definida).
 * @param b vetor de termos independentes.
 * @param M_inv matriz pré-condicionadora inversa (M⁻¹). Se NULL, usa a Identidade.
 * @param x vetor solução (deve ser pré-alocado e inicializado, geralmente com zeros).
 * @param n dimensão do sistema.
 * @param maxit número máximo de iterações.
 * @param epsilon tolerância para o critério de parada.
 * @param tempo_iter ponteiro para armazenar o tempo total das iterações.
 * @param norma_final ponteiro para armazenar a norma do erro na última iteração.
 * @return int O número de iterações executadas. Retorna -1 se não convergir.
 */
int solvePCG(real_t *A, real_t *b, real_t *M_inv, real_t *x, int n, long int maxit, 
    real_t epsilon, rtime_t *tempo_iter, real_t *norma_final);

#endif // __PCGC_H__