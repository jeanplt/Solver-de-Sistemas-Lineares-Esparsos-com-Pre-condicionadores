#!/bin/bash

# Entrada do script
if [ "$#" -ne 1 ]; then
    echo "Entrada errada"
    exit 1
fi

# Parâmetros para teste
label=""
n=0
k=0
maxit=0
epsilon="1e-8"

case "$1" in
    -l)
        label="leve"
        n=100
        k=5
        maxit=500
        ;;
    -m)
        label="media"
        n=1000
        k=11
        maxit=3000
        ;;
    -g)
        label="grande"
        n=2000
        k=15
        maxit=6000
        ;;
    -r)
        rm -f resultados_*.txt
        make purge
        exit 0
        ;;
    *)
        echo "Diretiva errada"
        exit 1
        ;;
esac

# Compilação do programa
if ! make; then
    echo "Erro na compilação"
    exit 1
fi

# Definição dos Omegas e Arquivo de Saída
omegas=(-1.0 0.0 1.0 1.5)
output_file="resultados_${label}.txt"

# Cria o arquivo de saída
echo "Valores de teste $label: n=$n k=$k omega(variável) maxit=$maxit epsilon=$epsilon" > "$output_file"
echo "" >> "$output_file"

# Executa os Testes alterando apeans o valor de Omega
for w in "${omegas[@]}"; do
    input="$n $k $w $maxit $epsilon"
    echo "Omega = $w:" >> "$output_file"
    echo "$input" | ./cgSolver >> "$output_file"
    echo "" >> "$output_file"
done

echo "Testes concluídos e salvos em $output_file"
make purge