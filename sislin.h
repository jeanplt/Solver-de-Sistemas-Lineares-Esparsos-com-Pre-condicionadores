#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

/**
 * Cria matriz 'A' k-diagonal e Termos independentes B
 * @param n dimensão da matriz A (tamanho n x n) e do vetor B (tamanho n)
 * @param k numero de diagonais da matriz A
 * @param A ponteiro para a matriz A (alocada dentro da função)
 * @param B ponteiro para o vetor B (alocado dentro da função)
 */
void criaKDiagonal(int n, int k, real_t **A, real_t **B);

/**
 * Gera a matriz simetrica positiva de A
 * @param A matriz k-diagonal original (tamanho n x n)
 * @param b vetor de termos independentes original (tamanho n)
 * @param n dimensão da matriz A e do vetor b
 * @param k numero de diagonais da matriz A
 * @param ASP ponteiro para a matriz simétrica positiva (alocada dentro da função)
 * @param bsp vetor de termos independentes da matriz simétrica positiva
 * @param tempo ponteiro para armazenar o tempo gasto na geração
 */
void genSimetricaPositiva(real_t *A, real_t *b, int n, int k, real_t **ASP, real_t *bsp, rtime_t *tempo);

/**
 * Devolve as matrizes D, L, U a partir da matriz A
 * @param A matriz k-diagonal original (tamanho n x n)
 * @param n dimensão da matriz A
 * @param k numero de diagonais da matriz A
 * @param D ponteiro para a matriz diagonal D (alocada dentro da função)
 * @param L ponteiro para a matriz estritamente inferior L (alocada dentro da função)
 * @param U ponteiro para a matriz estritamente superior U (alocada dentro da função)
 * @param tempo ponteiro para armazenar o tempo gasto na geração
 */
void geraDLU (real_t *A, int n, int k, real_t **D, real_t **L, real_t **U, rtime_t *tempo);

/**
 * Devolve a matriz M⁻¹
 * @param D matriz diagonal D (tamanho n x n)
 * @param L matriz estritamente inferior L (tamanho n x n)
 * @param U matriz estritamente superior U (tamanho n x n)
 * @param w pré-condicionador a ser utilizado (omega)
 * @param n dimensão das matrizes D, L, U
 * @param k numero de diagonais da matriz A original
 * @param M ponteiro para a matriz pré-condicionadora inversa M⁻¹ (alocada dentro da função)
 * @param tempo ponteiro para armazenar o tempo gasto na geração
 */
void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k, real_t **M, rtime_t *tempo);

/**
 * Calcula o resíduo || b - A*x ||l2
 * @param A matriz tamanho n x n (pode ser banda com zeros fora da banda)
 * @param b vetor RHS (tamanho n)
 * @param X solução candidata (tamanho n)
 * @param n dimensão
 * @param k número de diagonais (ímpar). Se k<=0 assume matriz densa (varre toda a linha).
 * @param tempo ponteiro para armazenar o tempo gasto na geração
 * @return norma L2 do resíduo (real_t)
 */
real_t calcResiduoSL (real_t *A, real_t *b, real_t *X, int n, int k, rtime_t *tempo);

#endif // __SISLIN_H__

