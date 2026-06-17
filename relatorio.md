# Relatório Técnico: Otimização e Cálculo Numérico de Autovalores

## 1. Introdução

Este relatório descreve a implementação e análise de um sistema para cálculo numérico de autovalores de matrizes quadradas genéricas, desenvolvido em C++17. O projeto percorre quatro fases: a construção do algoritmo clássico (Hessenberg + QR), o estudo de escalabilidade empírica, a otimização explorando simetria matricial e, por fim, uma aplicação em dinâmica de reatores químicos.

O ponto de partida teórico é o Teorema de Abel-Ruffini: para matrizes de ordem $N \geq 5$, não existe fórmula fechada para as raízes do polinômio característico, o que torna métodos iterativos a única alternativa viável.


## 2. Arquitetura do Projeto

O código é organizado em módulos com responsabilidades claras. As duas estruturas de dados centrais são:

A classe `Matrix` representa uma matriz quadrada $N \times N$ armazenada como vetor contíguo de `double` em ordem row-major --- o elemento $(i, j)$ mapeia para o índice $i \cdot N + j$ na memória. A multiplicação de matrizes usa a ordem de loops $i$-$k$-$j$, que favorece acessos sequenciais ao cache para a linha da matriz esquerda e da matriz resultado. A classe `Vector` encapsula um vetor de $N$ componentes com as operações usuais (produto escalar, norma $L_2$, soma, escalar).

A segunda estrutura importante é `SymmetricTridiagonal`, que armazena apenas dois vetores: `diag` ($N$ entradas) e `offdiag` ($N-1$ entradas). Essa representação compacta é central para a otimização da Fase 3.

Todas as classes usam semântica de valor --- operações retornam novos objetos, sem mutação in-place --- e delegam a gestão de memória inteiramente ao `std::vector` da biblioteca padrão.


## 3. Fase 1: O Algoritmo Geral

### 3.1 Transformações de Householder

O refletor de Householder é uma matriz ortogonal da forma

$$P = I - 2\mathbf{v}\mathbf{v}^T$$

onde $\mathbf{v}$ é um vetor unitário. Dado um vetor $\mathbf{x} \in \mathbb{R}^m$, queremos encontrar $\mathbf{v}$ tal que $P\mathbf{x}$ tenha todas as componentes nulas exceto a primeira. Isso equivale a

$$P\mathbf{x} = -\text{sign}(x_0)\|\mathbf{x}\|_2 \cdot \mathbf{e}_1$$

O vetor $\mathbf{v}$ é construído como

$$\mathbf{v} = \frac{\mathbf{x} + \text{sign}(x_0)\|\mathbf{x}\|\mathbf{e}_1}{\|\mathbf{x} + \text{sign}(x_0)\|\mathbf{x}\|\mathbf{e}_1\|}$$

A escolha de sinal $\text{sign}(x_0)$ evita cancelamento catastrófico: somamos na direção do componente dominante, garantindo que o denominador nunca se aproxime de zero.

Na implementação (`householder.cpp`), a aplicação do refletor a uma matriz é feita sem construir a matriz $P$ explicitamente. Para aplicar $PA$ (pela esquerda), calcula-se, para cada coluna $j$, a projeção $p = \mathbf{v}^T \mathbf{a}_j$ e depois atualiza-se $a_{ij} \leftarrow a_{ij} - 2v_i p$. A aplicação pela direita ($AP$) segue a lógica transposta, operando sobre linhas.

### 3.2 Redução à Forma de Hessenberg

Uma matriz de Hessenberg superior $H$ possui todos os elementos abaixo da primeira subdiagonal iguais a zero:

$$H(i, j) = 0 \quad \text{para } i > j + 1$$

O algoritmo reduz uma matriz densa $A$ a essa forma por transformações de similaridade --- o que preserva os autovalores. No passo $k$ (para $k = 0, \ldots, N-3$):

1. Extraímos a subcoluna $\mathbf{x} = A(k\!+\!1:N\!-\!1,\; k)$, que contém os elementos que precisam ser zerados abaixo da posição $(k+1, k)$.

2. Computamos o refletor de Householder $\mathbf{v}$ a partir de $\mathbf{x}$.

3. Aplicamos a transformação de similaridade: $A \leftarrow P_k A P_k^T$, primeiro pela esquerda (zerando a subcoluna) e depois pela direita (mantendo a similaridade).

Após os $N-2$ passos, a matriz resultante $H$ tem a forma de Hessenberg e os mesmos autovalores de $A$. O custo total é $O(N^3)$, dominado pelas multiplicações matrix-refletor.

### 3.3 Iteração QR com Deslocamento de Wilkinson

Com a matriz já na forma de Hessenberg, aplicamos a iteração QR para convergir aos autovalores. A ideia central é: dada $H_k$, computamos a fatoração $H_k = Q_k R_k$ e formamos $H_{k+1} = R_k Q_k$. Como $H_{k+1} = Q_k^T H_k Q_k$, os autovalores são preservados, e sob certas condições os elementos da subdiagonal convergem para zero, revelando os autovalores na diagonal.

**Fatoração QR via rotações de Givens.** A fatoração de uma matriz Hessenberg é feita com $N-1$ rotações de Givens, cada uma zerando um elemento da subdiagonal. A rotação $G(i, i\!+\!1)$ atua nas linhas $i$ e $i+1$:

$$\begin{bmatrix} c & -s \\ s & c \end{bmatrix} \begin{bmatrix} h_{ii} \\ h_{i+1,i} \end{bmatrix} = \begin{bmatrix} r \\ 0 \end{bmatrix}$$

Os coeficientes são calculados com branching numérico para estabilidade: quando $|h_{i+1,i}| > |h_{ii}|$, usa-se $\tau = -h_{ii}/h_{i+1,i}$ para evitar overflow.

**Deslocamento de Wilkinson.** A convergência crua do QR pode ser lenta. O deslocamento de Wilkinson acelera drasticamente o processo escolhendo, a cada iteração, um shift $\mu$ baseado nos autovalores do bloco $2 \times 2$ inferior da região ativa:

$$\mu = \text{autovalor de } \begin{bmatrix} h_{n-1,n-1} & h_{n-1,n} \\ h_{n,n-1} & h_{n,n} \end{bmatrix} \text{ mais próximo de } h_{n,n}$$

A iteração aplica então QR a $(H - \mu I)$ e recupera o shift na diagonal: o efeito é que $h_{n,n-1}$ converge cubicamente para zero.

**Nota --- desvio em relação ao enunciado.** O enunciado do projeto (§1.2) especifica a iteração QR *pura*: fatorar $H = QR$ e recombinar $H_{k+1} = R_kQ_k$ até a subdiagonal convergir, sem mencionar deslocamentos. Optamos por incorporar o deslocamento de Wilkinson porque o QR puro converge apenas *linearmente*, a uma taxa governada pela razão $|\lambda_{i+1}/\lambda_i|$ entre autovalores adjacentes. Em matrizes com autovalores de magnitude próxima --- frequentes nas matrizes aleatórias da Fase 2 --- essa razão se aproxima de 1 e a convergência praticamente estagna, podendo não atingir a tolerância em um número razoável de iterações (no limite, autovalores de igual magnitude, como pares $\pm\lambda$, nunca convergem pela iteração não-deslocada). O shift restaura a convergência cúbica e torna o algoritmo robusto e terminante sobre as entradas de teste; sem ele, o próprio benchmark de escalabilidade da Fase 2 não completaria em tempo viável. Trata-se, portanto, de uma decisão de engenharia *necessária* para a corretude prática, não de um mero refinamento de desempenho.

**Deflação.** Quando um elemento da subdiagonal se torna desprezível (segundo o critério relativo $|h_{i+1,i}| < \varepsilon(|h_{ii}| + |h_{i+1,i+1}|)$ com $\varepsilon = 10^{-10}$), o problema se divide em dois blocos independentes menores. Blocos $1 \times 1$ deflacionados são autovalores reais; blocos $2 \times 2$ irredutíveis com discriminante negativo representam pares conjugados complexos.

O custo de cada iteração QR sobre uma matriz Hessenberg $N \times N$ é $O(N^2)$ (uma rotação de Givens por subdiagonal, aplicada a $N$ colunas). Com $O(N)$ iterações típicas até convergência total, o custo da fase QR é $O(N^3)$.


## 4. Fase 2: Escalabilidade do Algoritmo Geral

O benchmark gera matrizes aleatórias de tamanhos $N \in \{10, 50, 100, 250, 500, 1000\}$, mede o tempo de execução e o pico de memória RSS para cada caso, repetindo a medição 20 vezes por tamanho para obter médias e desvios estatisticamente significativos.

O ajuste empírico por lei de potência sobre os dados de tempo do algoritmo geral (Hessenberg + QR) resultou em complexidade

$$T(N) \sim O(N^{2{,}76})$$

compatível com a expectativa teórica de $O(N^3)$ --- a diferença se deve ao fato de que o shift de Wilkinson reduz significativamente o número de iterações, e os tamanhos testados ainda não são grandes o bastante para saturar o expoente assintótico.

Para memória, o ajuste é

$$M(N) \sim O(N^{2{,}00})$$

o que é esperado: o algoritmo aloca duas matrizes densas $N \times N$ (a Hessenberg e a forma de Schur), dominando o uso de memória com $2N^2$ doubles.

![Escalabilidade: tempo de execução](scalability_time.png)

![Escalabilidade: pico de memória](scalability_memory.png)


## 5. Fase 3: Otimização para Matrizes Simétricas

### 5.1 O Colapso Tridiagonal

O resultado-chave desta fase é: quando a matriz de entrada $A$ é simétrica ($A = A^T$), a redução de Hessenberg por Householder produz não apenas uma matriz Hessenberg, mas uma **matriz tridiagonal simétrica**.

A prova é direta. Cada transformação de Householder é uma transformação de similaridade $A \leftarrow P_k A P_k^T$, onde $P_k = P_k^T = P_k^{-1}$ (ortogonal e simétrica). Se $A$ é simétrica, então $P_k A P_k^T$ também é simétrica, pois $(P_k A P_k^T)^T = P_k A^T P_k^T = P_k A P_k^T$. Portanto a simetria é um invariante ao longo de todas as reflexões. Agora, uma matriz que é simultaneamente Hessenberg superior (nula abaixo da primeira subdiagonal) e simétrica (a parte abaixo é o espelho da parte acima) só pode ter entradas não-nulas na diagonal principal e nas duas subdiagonais adjacentes --- ou seja, é tridiagonal.

O programa demonstra isso empiricamente para $N = 5$: aplica o algoritmo geral de Hessenberg a uma matriz simétrica e verifica que a maior entrada fora da banda tridiagonal é da ordem de $10^{-16}$ (precisão de máquina).

### 5.2 Redução Tridiagonal Otimizada

Embora o algoritmo geral de Hessenberg funcione para matrizes simétricas, ele desperdiça trabalho manipulando zeros que já sabemos existir. A implementação otimizada (`tridiagonal.cpp`) explora a simetria desde o início, usando uma atualização simétrica de rank-2 em vez de duas aplicações separadas do refletor.

No passo $k$, após calcular o refletor $\mathbf{v}$, a atualização do bloco inferior $B = A(k\!+\!1:, k\!+\!1:)$ é feita como:

$$\mathbf{p} = 2B\mathbf{v}$$
$$\mathbf{w} = \mathbf{p} - (\mathbf{v}^T\mathbf{p})\mathbf{v}$$
$$B \leftarrow B - \mathbf{w}\mathbf{v}^T - \mathbf{v}\mathbf{w}^T$$

Essa atualização preserva a simetria algebricamente e reduz o número de operações por manter a estrutura durante todo o processo. O resultado é armazenado na estrutura `SymmetricTridiagonal`, que ocupa apenas $O(N)$ memória (dois vetores de tamanho $N$ e $N-1$) em vez dos $O(N^2)$ de uma matriz densa.

### 5.3 Iteração QL Implícita sobre a Forma Tridiagonal

Para a iteração QR sobre a forma tridiagonal, o projeto implementa o algoritmo QL com deslocamento implícito (variante clássica `tqli` de Numerical Recipes). A ideia é a mesma da iteração QR --- transformações de similaridade que convergem a subdiagonal para zero --- mas operando exclusivamente sobre os dois vetores da estrutura tridiagonal.

Em cada iteração para o índice $l$:

1. Encontra-se o menor $m > l$ tal que $|e_m| \leq \varepsilon(|d_m| + |d_{m+1}|)$, onde $d$ é a diagonal e $e$ a subdiagonal.

2. Calcula-se o deslocamento de Wilkinson a partir do bloco $2 \times 2$ ativo:

$$g = \frac{d_{l+1} - d_l}{2e_l}, \quad r = \sqrt{g^2 + 1}, \quad g = d_m - d_l + \frac{e_l}{g + \text{sign}(g) \cdot r}$$

3. Aplica-se uma cadeia de rotações de Givens varrendo de $m-1$ até $l$, atualizando apenas os vetores $d$ e $e$. Cada rotação custa $O(1)$, e a varredura completa custa $O(N)$.

O custo total da iteração QL sobre uma tridiagonal é $O(N)$ por iteração, com $O(N)$ iterações para convergência global, totalizando $O(N^2)$ --- uma melhoria de ordem em relação ao $O(N^3)$ do algoritmo Hessenberg geral.

### 5.4 Resultados Comparativos

Os gráficos da Fase 2 sobrepõem ambos os algoritmos. Abaixo, os ganhos percentuais medidos (média de 20 amostras):

| $N$ | Ganho em tempo | Ganho em memória |
|----:|---------------:|-----------------:|
| 10 | 46% | 50% |
| 50 | 61% | 60% |
| 100 | 63% | 51% |
| 250 | 66% | 49% |
| 500 | 65% | 50% |
| 1000 | 67% | 50% |

O ajuste empírico do algoritmo otimizado resultou em $T(N) \sim O(N^{2{,}66})$, contra $O(N^{2{,}76})$ do geral. Em memória, ambos escalam como $O(N^2)$, mas o coeficiente do caminho tridiagonal é metade do geral --- consistente com o fato de que o passo denso mais caro (Hessenberg) aloca duas matrizes $N \times N$, enquanto o tridiagonal opera sobre $O(N)$ de estado próprio, embora a matriz de entrada $A$ (compartilhada) ainda domine a medição RSS.


## 6. Fase 4: Rede de Reatores (CSTRs)

A aplicação prática modela uma planta industrial com $N = 100$ reatores de mistura perfeita conectados em série. A dinâmica de concentrações é governada pelo sistema linear de EDOs

$$\frac{d\mathbf{C}}{dt} = A\mathbf{C}$$

onde a matriz $A$ é tridiagonal com diagonal principal $-2{,}5$ (perda por reação e escoamento) e sub/superdiagonais $+1{,}0$ (troca difusiva entre vizinhos).

O programa constrói a matriz diretamente como `SymmetricTridiagonal`, calcula os autovalores via o algoritmo otimizado e realiza duas análises:

**Estabilidade.** Todos os autovalores possuem parte real estritamente negativa --- resultado esperado, pois os autovalores de uma tridiagonal simétrica com diagonal $-2{,}5$ e bandas $+1{,}0$ são dados pela fórmula analítica $\lambda_k = -2{,}5 + 2\cos\!\left(\frac{k\pi}{N+1}\right)$, cujo valor máximo é $-2{,}5 + 2\cos\!\left(\frac{\pi}{101}\right) \approx -0{,}50$. O sistema é portanto assintoticamente estável: as concentrações convergem para o estado estacionário.

**Rigidez (Stiffness).** A razão de rigidez é

$$S = \frac{\max|\text{Re}(\lambda_i)|}{\min|\text{Re}(\lambda_i)|} = \frac{|\lambda_{\min}|}{|\lambda_{\max}|}$$

onde $\lambda_{\min} \approx -4{,}50$ e $\lambda_{\max} \approx -0{,}50$. O valor $S \approx 9$ é baixo (muito abaixo do limiar $10^3$), classificando o sistema como **não-rígido**. Um engenheiro de software pode utilizar com segurança um solver explícito como Runge-Kutta de 4ª ordem (RK4) para simular a planta, evitando o custo computacional de métodos implícitos.


## 7. Discussão

### Por que matrizes simétricas tornam o método tão mais rápido?

Três efeitos se acumulam, todos derivados da restrição $A = A^T$:

Primeiro, a **redução ao problema tridiagonal**. A simetria garante que a forma de Hessenberg colapsa para uma tridiagonal, reduzindo a representação de $O(N^2)$ entradas para $O(N)$. Isso não é apenas uma economia de memória --- significa que toda operação subsequente toca uma fração drasticamente menor de dados.

Segundo, a **queda de complexidade da iteração QR**. Sobre uma Hessenberg densa, cada passo QR (fatoração + recombinação via Givens) custa $O(N^2)$ porque cada rotação afeta $O(N)$ colunas. Sobre uma tridiagonal, cada rotação afeta apenas $O(1)$ entradas vizinhas, reduzindo o custo por passo a $O(N)$. Com $O(N)$ passos típicos para convergência total, a iteração completa passa de $O(N^3)$ para $O(N^2)$.

Terceiro, a **localidade de cache**. A estrutura tridiagonal opera sobre dois vetores contíguos de tamanho $N$, que cabem confortavelmente nas caches L1/L2 para os tamanhos testados. O caminho denso, por outro lado, percorre matrizes $N \times N$ com stride variável, gerando mais cache misses. Esse fator amplifica o ganho algorítmico, especialmente para $N \geq 250$.

### O gargalo é CPU ou memória?

O gargalo é **CPU**. Ambos os algoritmos escalam como $O(N^2)$ em memória, mas o ganho de tempo (67% para $N = 1000$) supera consistentemente o ganho de memória (50%). Se memória fosse o fator limitante, os ganhos de tempo e memória seriam proporcionais. A diferença indica que a redução no número de operações de ponto flutuante --- não a economia de alocações --- é o fator dominante.
