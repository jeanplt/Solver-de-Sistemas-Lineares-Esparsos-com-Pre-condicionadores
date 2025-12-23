# Solução de Sistemas Lineares Esparsos com Pré-condicionadores (Geral)

# Autoria do software

    Autor:      Jean Pablo Lopes Teixeira
    GRR:        20224382
 
# Apresentação geral

    Este programa modular segue os parâmetros especificados pelo enunciado do 
    trabalho, no qual temos a aplicação do método dos gradientes conjugados 
    para resolução de sistemas lineares esparsos com pré-condicionadores. 

# Componentes

    Dos arquivos que compõem este programa modular, temos as seguintes 
    nomenclaturas e suas funções:

    # utils.h e utils.c: Respectivamente a biblioteca de definições gerais
    como obtenção de valor absoluto, valores reais, inteiros, char/strings e
    tempo, bem como protótipos para função de marcação de tempo e contador 
    para likwid, e o arquivo com as funções da biblioteca implementadas.

    # sislin.h e sislin.c: Conjunto biblioteca e arquivo de implementação
    base para o desenvolvimento do programa, com todas as funções necessárias
    para gerar a matriz A com termos independentes em B, gerar uma matriz
    simétrica positiva, gerar a partir de uma matriz A outras três matrizes,
    sendo elas matriz triangular inferior, matriz triangular superior e matriz
    diagonal principal, função para gerar o pré-condicionador a partir da fórmula
    M = (D + ωL)D⁻¹ (D + ωU ), ω ∈ [0, 2], retornando M⁻¹, além de uma função
    para o cálculo do resíduo (||r||l2), onde r = b - Ax.

    # pcgc.h e pcgc.c: A biblioteca apresenta todos os protótipos das funções
    que serão usadas para calcular o produto escalar entre dois vetores, 
    atualizar os vetores de solução e resíduo a cada iteração, multiplicar
    a matriz por vetor e calcular a norma máxima para o critério de parada. 
    Além destas funções auxiliares, temos a função chave deste conjunto,
    responsável por resolver o sistema linear Ax=b usando o método dos 
    Gradientes Conjugados Pré-condicionados, com utilização de todas as
    funções supracitadas.

    # cgSolver.c: Arquivo principal do programa modular, possui as entradas
    solicitadas de dados, seguido da inicialização de todas as matrizes e
    vetores necessários para geração do sistema k-diagonal, transformação 
    para um sistema simétrico positivo, geração do pré-condicionador M⁻¹,
    resolução do sistema, cálculo do resíduo e impressão dos resultados, com
    o devido registro de tempo para cada atividade e liberação de memória do
    programa principal.

    # Makefile: Ferramenta de compilação do programa modular, com toda a base
    de compilação e regras base para atender aos requisitos do trabalho.

    # LEIAME: Este arquivo, do qual está redigida toda a documentação 
    introdutória do programa modular.

    # Shell.sh : Script que não agrega ao programa modular, porém é útil para 
    realizar testes em escala com os quatro tipos de omegas disponíveis e salvar
    seus resultados em um arquivo .txt. Para rodá-lo basta realizar a compilação 
    básica de shell com chmod +x Shell.sh e executá-lo com uma das quatros 
    diretivas abaixo, no seguinte formato:
        
        ./Shell.sh -l: Realiza o teste leve, com sistema linear pequeno.
        ./Shell.sh -m: Realiza o teste médio, com sistema linear intermediário.
        ./Shell.sh -g: Realiza o teste grande, com sistema linear maior.
        ./Shell.sh -r: Remove os arquivos .txt gerados dos testes acima.
    
    Observação: O teste grande tende a demorar alguns segundos, que podem se 
    prolongar para minutos se n ou k forem aumentados ainda mais.

# Estrutura de Dados

    Da estrutura de dados presente neste programa modular, no geral temos 
    estruturas mais simples que possuem uma clara separação de 
    responsabilidades, sendo esta separação da seguinte forma:

    # Tipos de dados base (typedef): Definidos em utils.h e utilizados ao
    longo de todo o programa, principalmente para cálculo de valores com
    ponto flutuante pelo uso de real_t, valores grandes de inteiros com
    size_t e também intervalos de tempo com rtime_t.

    # Arrays: Alocados dinamicamente com tamanho size_t, representam os
    vetores de termos independentes das matrizes presentes em sislin.c e
    cgSolver.c, o vetor solução em cgSolver.c e todos os vetores auxiliares
    para cálculo do sistema linear pela função chave de pcgc.c.

    # Matrizes: Também alocadas dinamicamente com tamanho size_t², e expressas
    de forma linearizada para facilitar o acesso de dados. Estão presentes na
    inicialização do programa em cgSolver, na geração de matriz em sislin.c e
    em todas as operações que envolvam a resolução do sistema linear em 
    cgSolver.c, além de serem desalocadas tanto no cenário de erro quanto de 
    sucesso da rodagem do programa pela main.

# Modularização

    A modularidade deste programa é obtida através dos seguintes arquivos que 
    estão interligados pelas suas bibliotecas de protótipos:

    # Módulo principal: Em cgSolver.c, manipulamos quase toda a memória 
    principal do programa por meio de matrizes e vetores, dos quais declaramos 
    seus ponteiros para serem utilizados pelas funções de sislin.c e pcgc.c, e 
    após o preenchimento destes ponteiros, escrevemos os resultados necessários 
    e liberamos toda a memória utilizada.

    # Módulo de sistema linear: Em sislin.c, são feitas as alocações e 
    preparações dos parâmetros recebidos, sendo estes as matrizes necessárias 
    para criação de uma matriz k-diagonal, simétrica positiva e também 
    pré-condicionadores, além da geração dos vetores de termos independentes que 
    estarão todos preparados para utilização dentro do módulo principal e módulo 
    de solução após cada chamada.

    # Módulo de solução: Em pcgc.c, realizamos todo o cálculo e processamento do 
    sistema linear, com funções auxiliares que fazem todos os cálculos a partir 
    de valores inicializados pelas funções de sislin.c, e também vetores 
    auxiliares alocados na função chave para manipular os resultados obtidos e 
    trazer à tona uma solução útil para o sistema ou, no pior caso, um erro de 
    não convergência. Ao final das operações, retornamos ao módulo principal para 
    encerrar o programa. 
