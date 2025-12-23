#include "utils.h"
#include "sislin.h"

static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k );
static inline real_t generateRandomB( unsigned int k );

/**
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * @param i,j coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return ( (i==j) ? (real_t)(k<<1) : 1.0 )  * (real_t)random() * invRandMax;
}

/**
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomB( unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return (real_t)(k<<2) * (real_t)random() * invRandMax;
}

void criaKDiagonal(int n, int k, real_t **A, real_t **B) 
{
  /* Inicializações de inteiros i e j para loops, matriz A e vetor B */ 
  int i, j;
  size_t nL = (size_t)n; 
  real_t *a = calloc(nL*nL, sizeof(real_t));
  real_t *b = calloc(nL, sizeof(real_t));

  if (!a || !b) {
    fprintf(stderr, "Erro na alocação de memória na criaKDiagonal\n");
    return;
  }

  int half_k = k / 2;

  /* Laço para gerar a matriz A e vetor B com gerador aleatório */
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
        /* Verifica se está dentro da banda para gerar valor aleatório, se não força zero */
        if (ABS(i - j) <= half_k) {
            a[i * n + j] = generateRandomA(i, j, k);
        } else {
            a[i * n + j] = 0.0;
        }
    }
    b[i] = generateRandomB(k);
  }

  /* Carrega a matriz a e vetor b nos parâmetros */
  *A = a;
  *B = b;
} 

void genSimetricaPositiva(real_t *A, real_t *b, int n, int k, real_t **ASP, real_t *bsp, rtime_t *tempo)
{
  *tempo = timestamp();
  
  /* Inicializações de inteiros i e j para loops, half para largura de banda e asp para matriz */
  int i, j, p;
  int half = k / 2;
  size_t nL = (size_t)n; 

  real_t *asp = calloc(nL * nL, sizeof(real_t));
  if(!asp) {
    fprintf(stderr, "Erro na alocação de memória na genSimetricaPositiva\n");
    *tempo = timestamp() - *tempo;
    return;
  }
  
  /* Laço para calcular asp = A * A^T*/
  for (i = 0; i < n; i++) {
    size_t row_i_idx = (size_t)i * nL;
    
    /* Calcula apenas j >= i devido à simetria e banda aumentada de asp */
    int jmax = (i + 2 * half > n - 1) ? n - 1 : i + 2 * half;

    for (j = i; j <= jmax; j++) {
      /* Para acesso direto na matriz N x N */
      size_t row_j_idx = (size_t)j * nL;
      real_t sum = 0.0;

      /* Calcula banda da linha i: [i - half, i + half] */
      int start_i = (i - half < 0) ? 0 : i - half;
      int end_i   = (i + half > n - 1) ? n - 1 : i + half;

      /* Calcula banda da linha j: [j - half, j + half] */
      int start_j = (j - half < 0) ? 0 : j - half;
      int end_j   = (j + half > n - 1) ? n - 1 : j + half;

      /* Calcula a interseção */
      int p_start = (start_i > start_j) ? start_i : start_j;
      int p_end   = (end_i < end_j)     ? end_i   : end_j;

      /* Verifica intervalo para calcular o produto escalar */
      if (p_start <= p_end) {
          for (p = p_start; p <= p_end; p++) {
            sum += A[row_i_idx + (size_t)p] * A[row_j_idx + (size_t)p];
          }
      }
      
      /* Verifica positividade na diagonal principal */
      if (i == j)
        sum += (real_t)k;
      
      /* Armazena o resultado em asp e verifica simetria */
      asp[row_i_idx + (size_t)j] = sum;
      if (i != j) 
          asp[row_j_idx + (size_t)i] = sum;
    }
  }

  /* Laço para calcular bsp = ASP * b */
  for (i = 0; i < n; ++i) {
    long double sum = 0.0L;
    size_t row_idx = (size_t)i * nL;
    
    /* Define limete da banda na matriz asp */
    int jmin = (i - 2 * half < 0) ? 0 : i - 2 * half;
    int jmax = (i + 2 * half > n - 1) ? n - 1 : i + 2 * half;

    for (j = jmin; j <= jmax; ++j) 
      sum += (long double)asp[row_idx + (size_t)j] * (long double)b[j];
    
    
    bsp[i] = (real_t)sum;
  }
  
  /* Carrega a matriz asp no parâmetro */
  *ASP = asp;
  *tempo = timestamp() - *tempo;
}

void geraDLU (real_t *A, int n, int k, real_t **D, real_t **L, real_t **U, rtime_t *tempo)
{
  *tempo = timestamp();

  /* Inicializações de inteiros i e j para loops, matriz D, L e U para A */
  int i, j;
  size_t nL = (size_t)n; // Para evitar overflow na alocação com multilplicação
  real_t *d = calloc(nL*nL, sizeof(real_t));
  real_t *l = calloc(nL*nL, sizeof(real_t));
  real_t *u = calloc(nL*nL, sizeof(real_t));

  if (!d || !l || !u) {
      fprintf(stderr, "Erro na alocação de memória na geraDLU\n");
      *tempo = timestamp() - *tempo;
      return;
  }

  /* Calcula half = k/2 e limita o valor para [0, n-1] */
  int half = (k > 0) ? (k / 2) : (n - 1);
  if (half > n-1) 
    half = n-1;

  /* Laço para percorrer A e preencher D, L, U */
  for (i = 0; i < n; ++i) {
      size_t row_i = (size_t)i * nL;

      /* Calcula os limites de banda para linha */
      int jmin = i - half;
      if (jmin < 0) 
        jmin = 0;
      
      int jmax = i + half;
      if (jmax > n - 1) 
        jmax = n - 1;

      /* Laço para preencher d, l, u conforme posição */
      for (j = jmin; j <= jmax; ++j) {
          real_t aij = A[row_i + (size_t)j];
          if (i == j) {
              d[row_i + (size_t)j] = aij;   // diagonal
          } else if (i > j) {
              l[row_i + (size_t)j] = aij;   // estritamente inferior
          } else { // i < j 
              u[row_i + (size_t)j] = aij;   // estritamente superior
          }
      }
  }

  /* Carrega as matrizes d, l e u nos parâmetros */
  *D = d;
  *L = l;
  *U = u;

  *tempo = timestamp() - *tempo;
}

void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k, real_t **M, rtime_t *tempo)
{
  if (tempo) *tempo = timestamp();

  /* Inicializações de inteiros para loops, matriz M e vetores z, y e x */
  int i, j, p, q, row;
  size_t nL = (size_t)n; // Para evitar overflow na alocação com multilplicação
  real_t *Minv = calloc(nL * nL, sizeof(real_t));
  real_t *z  = malloc(nL * sizeof(real_t));
  real_t *y  = malloc(nL * sizeof(real_t));
  real_t *x  = malloc(nL * sizeof(real_t));

  if (!Minv || !z || !y || !x) {
      fprintf(stderr, "Erro na alocação de memória na geraPreCond\n");
      if (tempo) *tempo = timestamp() - *tempo;
      free(Minv);
      free(z); free(y); free(x);
      return; 
  }

  /* Calcula half = k/2 e limita o valor para [0, n-1] */
  int half = (k > 0) ? (k / 2) : (n - 1);
  if (half > n-1) 
    half = n-1;

  // Pré-cálculo do tTmanho do vetor z para o memset
  size_t z_bytes = nL * sizeof(real_t);

  /* Laço para calcular cada coluna j de M⁻¹ */
  for (j = 0; j < n; ++j) {
    
    /* Zera o vetor z a cada iteração para substituição direta */
    memset(z, 0, z_bytes);

    /* Laço para substituição direta (D + w L) z = e_j */
    for (i = 0; i < n; ++i) {
        /* Acumulador s como long double para maior precisão */
        long double s = 0.0L;
        int jmin = i - half;
        if (jmin < 0) 
          jmin = 0;

        /* Pré-cálculo do início da linha i */
        size_t row_i = (size_t)i * nL;
        
        /* Laço para iterar dentro da banda */
        for (p = jmin; p < i; ++p) {
            real_t Lip = L[row_i + (size_t)p];
            s += (long double)Lip * (long double)z[p];
        }

        /* Para cálculo de zi e acesso a diagonal D */
        long double bi = (i == j) ? 1.0L : 0.0L;
        real_t Di = D[row_i + (size_t)i]; 

        if (fabsl((long double)Di) < 1e-300L) {
            fprintf(stderr, "1⁰ Divisao por zero na geraPreCond (i=%d)\n", i);
            free(Minv);
            free(z); free(y); free(x);
            if (tempo) *tempo = timestamp() - *tempo;
            return; 
        }
        
        /* Cálculo de z[i] */
        long double zi = (bi - (long double)w * s) / (long double)Di;
        z[i] = (real_t)zi;
        
        /* Cálculo de y[i] */
        real_t Dii = D[row_i + (size_t)i];
        y[i] = Dii * z[i];
    }

    /* Laço para substituição reversa (D + w U) x = y */
    for (i = n - 1; i >= 0; --i) {
        /* Acumulador s como long double para maior precisão */
        long double s = 0.0L;
        int jmax = i + half;
        if (jmax > n-1) 
          jmax = n-1;
        
        /* Pré-cálculo do início da linha i */
        size_t row_i = (size_t)i * nL;
        
        /* Laço para iterar dentro da banda */
        for (q = i + 1; q <= jmax; ++q) {
            real_t Uiq = U[row_i + (size_t)q];
            s += (long double)Uiq * (long double)x[q];
        }
        
        real_t Di = D[row_i + (size_t)i];

        if (fabsl((long double)Di) < 1e-300L) {
            fprintf(stderr, "2⁰ Divisao por zero na geraPreCond (i=%d)\n", i);
            free(Minv);
            free(z); free(y); free(x);
            if (tempo) *tempo = timestamp() - *tempo;
            return;
        }

        // Calcula xi e armazenamento em x[i]
        long double xi = ((long double)y[i] - (long double)w * s) / (long double)Di;
        x[i] = (real_t)xi;
    }

    /* Laço para armazenar coluna x (vetor solução) na matriz Minv */
    for (row = 0; row < n; ++row) {
        Minv[(size_t)row * nL + (size_t)j] = x[row];   
    }
  }

  /* Carrega a Matriz M e libera os vetores utilizados */
  *M = Minv;
  free(z); free(y); free(x);

  if (tempo) *tempo = timestamp() - *tempo;
}

real_t calcResiduoSL (real_t *A, real_t *b, real_t *X, int n, int k, rtime_t *tempo)
{
  *tempo = timestamp();

  /* Inicializações de inteiros i e j para loops, half e sumsq para cálculo de banda */ 
  int i, j;
  int half = k / 2;
  long double sumsq = 0.0L;

  /* Laço para calcular o resíduo || b - A*x ||l2 */
  for (i = 0; i < n; ++i) {
    /* Variáveis de intervalo de colunas onde os valores da matriz A
    são geralmente diferentes de zero */
    int jmin = i - half; 
    if (jmin < 0) 
      jmin = 0;
    
    int jmax = i + half; 
    if (jmax > n - 1) 
      jmax = n - 1;
    
    long double Ax_i = 0.0L;

    /* Variável e laço para multiplicação só dentro da banda */
    size_t row = (size_t)i * (size_t)n;
    for (j = jmin; j <= jmax; ++j) 
      Ax_i += (long double) A[row + (size_t)j] * (long double) X[j];

    long double ri = (long double) b[i] - Ax_i;
    sumsq += ri * ri;
  }

  /* Tira a raíz do resultado para resposta */
  real_t res = (real_t) sqrt((double) sumsq);

  *tempo = timestamp() - *tempo;
  return res;
}