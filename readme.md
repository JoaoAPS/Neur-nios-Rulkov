Documentação: Rulkov Lib Versão 3.0
===================================

*28/01/19*

Esta biblioteca implementa redes em que cada elemento é um neurônio de Rulkov, e cada conexão representa
uma sinapse entre dois neurônios.

Para mais informações sobre neurônios de Rulkov, ver [link](https://en.wikipedia.org/wiki/Rulkov_map).

Para instruções de como usar a biblioteca abra com seu navegador o arquivo "manual.html" na pasta "docs"

# Dinâmica Local

Os neurônios de Rulkov são definidos por duas variáveis: `x`, que representa o potencial de membrana, 
e `y`, uma variável auxiliar de evolução lenta.
A evolução temporal de `x`  apresenta bursts, vários spikes em sequência.

O modelo de Rulkov é um mapa, isto é, apresenta tempo discreto.


# Acoplamento

*Conexões* se refere a como cada elemento da rede está conectado com os outros.
Suas propriedades são definidas pela *estrutura* (quais elementos estão conectados com quais),
e pela forma com que elementos conectados se influenciam (aqui chamado *sinapses*).

## Estrutura

A estrutura dita quais neurônios estão conectados com quais.
Há uma divisão logística entre dois grupos de estruturas: global e não-globais.

Redes globais são aquelas em que todos elemetos estão conectados com todos os outros.
Para utilizar uma rede global, você deve criar um objeto `RedeGlobal`, presente em *"include/Global.h"*.

Já redes não-globais têm um conectoma arbitrário, sendo ele definido por uma matriz de adjacência `A_ij`, que é 1 se há conexão entre o neurônio $j$ e o $i$, e é 0 se não houver.
Mas para economia de memória, o programa não usa uma matrix de adjacência, mas sim um vetor de adjacência, que é um vetor
de inteiros, sendo que cada inteiro representa uma conexão.
Essa representação é dada por `adjVetElem = i * N + j`, onde `ì` e `j` são a linha e coluna, respectivamente, 
de um elemento não-nulo da matrix de adjacência.

Para utlilizar uma rede não-global, você deve criar um objecto `Rede`, presente em *"include/Rede.h"*, e passar um arquivo
com o vetor de adjacência para ser lido, através de
```c++
rede.readAdjVet( |adjVet filePath| );
```

O arquivo do vetor de adjacência deve estar formatado com um número por linha, obrigatoriamente em ordem crescente, 
cada um representando uma conexão.

## Sinapses
Sinapses representam a forma como dois neurônios conectados influenciam um ao outro.

Está presente na biblioteca o código para dois tipos de sinapses: químicas e elétricas,
embora atualmente apenas as químicas estejam de fato sendo usadas na biblioteca.

Ambos tipos de sinapses têm suas forças controladas pelo *parâmetro de acoplamento* (denotado `eps`).
Se seu valor for zero, os neurônios não se influenciam em nada; e quanto maior for `eps`, maior será a
interação entre os neurônios.

Além disso, elas também são normalizadas pelo tamanho da rede, para que os mesmos valores `eps`s tenham efeitos
similares em redes diferentes.
Há várias formas de fazer a normalização.
A aqui usada consiste em dividir a corrente sináptica pelo grau médio de conexões, isto é, a média do número de vizinhos 
de cada neurônio.

### Sinapses Elétricas
Atualmente não são usadas, mas seu código ainda está presente na biblioteca.

Essas são as sinapses de neurônios que estão em contato direto com o outro.
Isso significa que o potencial de membrana de um influencia diretamente o outro.
A forma da corrente sinaptica nesse caso é trivial, apenas proporcional ao potencial do neurônio vizinho.
```c++
correnteSinaptica = (eps / conectividadeMedia) * x_vizinho
```

### Sinapses Químicas
Representam a interação que neurônios realizam por neurotransmissores.

A forma usada é a mesma que pelos outros alunos do Viana nos artigos de 2018 e 2019,
como [este](http://fisica.ufpr.br/viana/artigos/2018/Mugnaine2018_Article_DelayedFeedbackControlOfPhaseS.pdf).
Ela não depende do potencial de membrana do vizinho, apenas do potencial do neurônio afetado.


A implementação da sinapse no programa é feita da seguinte forma:
Ela só é acionada quando o neurônio présinaptico (`j`) está disparando.
Isso é implementado através de uma função degral no potencial de `j` com limiar `theta`.
A corrente injetada é proporcional a
```c++
correnteSinaptica = (eps / conectividadeMedia) * (constanteSinaptica - x_pos)
```
`constanteSinaptica` é uma constante que determina se a sinapse é excitatória ou inibitória.
Atualmente apenas sinapses excitatórias são usadas.
`x_pos` é o valor instantâneo do potencial do neurônio que está recebendo a corrente sináptica (pós-sináptico).

Os valores das constantes são `theta = -1` e `V = 1`, conforme o artigo de 2018 do Fabiano Ferrari.


# Plasticidade

Plasticidade se refere à possibilidade das sinapses terem suas forças alteradas durante a evolução da rede.

Existem alguns protocolos que informam como os pesos de sinapses mudam em uma rede de neurônios.
Aqui implementamos dois deles: a STDP (Spike-Time Dependent Plasticity) - que altera os pesos de sinapse
com base na diferença de tempo entre spikes de dois neurônios conectados - e a BTDP (Burst-Time Dependent Plasticity) - 
que altera os pesos de sinapse com base na diferença de tempo entre os inícios de burst de dois neurônios conectados.

Embora a STDP seja mais estudada e firmada, ela trabalha com spikes.
Então a BTDP talvez seja a melhor alternativa, mesmo sendo relativamente mais nova e menos consolidada.

## Plasticidade BTDP

A BTDP foi introduzida por Butts, como uma plasticidade que modela os pesos de sinapse com base na defasagem de 
bursts de neurônios conectados.
Neurônios que iniciam seus bursts em tempos próximos têm sua força de conexão aumentada, enquanto bursts que 
distam temporalmente muito geram depressão do pesos de sinapse.
Vale notar que a ordem em os bursts ocorrem não influencia a mudança.
A função de BTDP é linear e simétrica, tem um pico em 0 e decresce tanto para o lado positivo quanto para o negativo linearmente, 
até um certo *tempo de saturação*, a partir do qual ela vira constante.

A implementação computacional é feita da seguinte forma:

Quando um neurônio `i` inicia um burst, corremos por todas suas ligações (que entram e que saem), e com base na diferença de tempo entre este bst (burst start time) do `i` e o último bst de cada vizinho (conforme a regra BTDP descrita acima).
Com esse algoritmo, uma sinápse entre `i` e `j` é atualizada em duas situações: 1) Quando `i` dispara; e 2) Quando `j` dispara.

Entretanto, como a alteração de peso de sinapse ocorre devido à interação entre dois bursts, gostaríamos que houvesse uma atualização
a cada par de bursts `i` e `j`.
Para resolver isso, a regra BTDP implementada no código é levemente diferente da de Butts, tendo a mesma forma mas amplitudes diferentes.
A saber, se a curva BTDP de Butts tem amplitude de depressão `D`, amplitude de potencialização `P`, e tempo de stauração `T`,
então a curva BTDP usada na implemetação tem amplitude de depressão `d = D/2`, amplitude de potencialização `p = P + d`, 
e tempo de stauração `T`.
Com essas amplitudes, e usando o algoritmo descrito no parágrafo anterior, a regra BTDP de Butts é obtida a cada par de 
bursts `i` `j`.
Segue um exemplo na prática.

```
Considere os seguintes tempos de burst:
i: 0, 105
j: 100

i:    -|-------------------------------------------------------|-
j:    -----------------------------------------------------|-----
tempo: 0                                                  100 105

Em t=100, j dispara, e a regra é acionada.
Como a diferença entre o tempo atual e o último tempo de burst de `i` é grande (100-0 = 100), 
ocorrerá depressão máxima e a sinápse será mudada de -d=-D/2.
Em t=105, i dispara, e a regra é acionada.
Como a diferença entre o tempo atual e o último disparo de j é pequena (105-100 = 5), 
então ocorrerá potencialização aproximadamente máxima, e a sinápse será mudada de p=P+d.

Ou seja, no final da interação, a sinapse foi alterada de -d+p = -d + P + d = P.
Que condiz com o esperado de que dois bursts próximos tenham potencialização P.

Se o segundo burst de i ocorresse muito depois de t=100, então novamente ocorreria depressão de -d=-D/2, 
de modo que a mudança total seria -d-d = -2 D/2 = -D.
```

Além da regra de mudança BTDP, forçamos que o valor do peso de sinapse esteja entre limitado, atribuindo valores máximo e mínimo.
No programa, se algum peso de sinapse ultrapassa esses valores, ele é trazido à mão ao valor limitante que ele ultrapassou.
Os valores mínimo e máximo usados são 0 e 1 (lembrando que todas correntes de sinapse são multiplicadas por `eps`, que 
representa o peso de sinápse máximo, tal que a sinapse efetiva está limitada entre 0 e `eps`).

Os parâmetros default usados são: 

- Amplitude de Potencialização: `P = 0.008`. Escolhido arbritariamente para que a mudança ocorra lentamente, mas não tanto.
- Razão entre amplitudes de depressão e potencialização: `R = D / P = 0.4`. Razão aproximada da curva de Butts.
- Tempo de Saturação: `T_s = 59 steps`. Escolhido para manter a mesma proporção `IBI/T_s` e `burstDuration/T_s` que as encontradas
para os neurônios OFF de Gjorgjieva (que foram o tipo de neurônio dos artigos biológicos mais parecidos com o Rulkov, no
sentido de razões entre tempos característicos.)


## Plasticidade STDP

Aqui consideramos um protocolo STDP simples, em que a interação entre cada spike do pré com cada spike do pós dentro de 
um burst é levada em conta.
Se o neurônio pré `j` tem um spike logo antes de um spike do pós `i`, a força de sinapse de `j` para `i` é aumentada por
`A_p * exp( -(t_i - t_j) / tau_p )`.
Já se o neurônio pré `j` dispara depois do pós `i`, a força de sinapse de `j` para `i` é diminuida por 
`A_d * exp( -(t_j - t_i) / tau_d )`.
Aqui `t_i` é o tempo de início de burst do neurônio `i`; e `A_p`, `A_d`, `tau_p` e `tau_d` são parâmetros constantes.

A implementação computacional é feita da seguinte forma:

Quando um neurônio `i` tem um spike, corremos por todas as ligações que **saem** dele e diminuimos a força dessas ligações
com base na diferença de tempo entre este spike do `i` e os *eta* últimos spikes de cada vizinho (conforme descrito acima).
Isso significa que todas as ligações que saem de `i` serão diminuídas, mas o decaimento exponencial dessa mudança com a 
diferença de tempo dos spikes faz com que essa mudança seja muito pequena quando os tempos de disparo dos neurônios `i` e
`j` estão muito separados temporalmente.
Em seguida corremos por todas que **chegam** em `i` e aumentamos a força de todas elas, novamente com base no decaimento 
exponencial descrito acima.

Além disso, forçamos que o valor da força de sinapse esteja entre limitado, atribuindo valores máximo e mínimo.
Se alguma força de sinapse ultrapassa esses valores, ela é trazida à mão ao valor limitante que ela ultrapassou.
Os valores mínimo e máximo usados são 0 e 1 (lembrando que todas correntes de sinapse são multiplicadas por `eps`, 
tal que a sinapse efetiva está limitada entre 0 e `eps`)
Esse tipo de técnica poderia ser substituído por parâmetros `A_p` e `A_d` que dependem da força de sinapse atual de 
maneira que ela é naturalmente limitada.

Os parâmetros default usados são: 

- Amplitude de Pontencialização `A_p = 0.0001`
- Amplitude de Depressão `A_d = 0.00008`
- Tempo de decaimento característico da Potencialização `tau_p = 20 (~20ms)`
- Tempo de decaimento característico da Depressão `tau_d = 20`.

A razão 0.8 entre A_d e A_p foi tirada de 'Zhang LI (Nature 395 1998) A critical window for cooperation a ...'

