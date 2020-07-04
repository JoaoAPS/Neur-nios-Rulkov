Manual: Rulkov Lib Versão 3.0
=============================

*28/01/19*

Estruturas Suportadas e o Arquivos em que se encontram

- Rede Global -> Global.cpp
- Rede com Arquitetura segundo um vetor de Adjacência -> Rede.cpp
- Rede segundo vetor de Adjacência com plasticidade nos pesos de sinápses -> Plasticidade.cpp
- Sinal Externo Harmônico a ser adicionado em uma das Redes acima -> SinalHarmonico.cpp



# MANUAL

## Compilação da biblioteca
Para compilar a biblioteca execute o comando `make` no diretório da biblioteca.
É necessário ter o `g++` instalado.
A compilação gera um arquivo "librulkov.a", que deve ser linkado na compilação do seu programa.

## Na compilação do programa

#### - Compile o arquivo .a junto com o .cpp:
`$(CC) main.cpp /path/to/libs/librulkov.a -o exec`

O arquivo .a deve vir depois do .cpp na chamada da compilação

Opcional: Incluir flag '-I/path/to/lib/rulkov/' para incluir headers com <|file|.h> no arquivo .cpp



## No seu arquivo .cpp

#### - Inclua os Headers necessários:	
```c++
#include "/path/to/libs/rulkov/Rede.h"
#include "/path/to/libs/rulkov/Global.h"
```
Apenas um desses é necessário.


#### - Crie o objeto Rede correspondente à estrutura desejada:
```c++
Rede	   rede(NUM_NEURONS, NUM_STEP, TRANSIENTE, [EPS, BETA, SIGMA]);
RedeGlobal rede(NUM_NEURONS, NUM_STEP, TRANSIENTE, [EPS, BETA, SIGMA]);
RedePlasticidade rede(NUM_NEURONS, NUM_STEP, TRANSIENTE, SYNWEI_MAX, MODELO);
```

Os valores default são:
	EPS   = 0.0
	BETA  = 0.001
	SIGMA = 0.001


#### - (OPCIONAL) Defina a seed de números aleatórios:
```c++
rede.setSeed( |seed| );
```


#### - Inicie a rede:
```c++
rede.init([alpha], [x0], [y0]);
```

+ alpha pode ser um float para que todos neurônios da rede tenham o mesmo valor de alpha;
+ alpha pode ser omitido para uma distribuição aleatória.

+ x0 e y0 podem ser objetos std::vector<float>.
	* Se o número de elementos desses vectors for igual ao número de neurônios, cada neurônio da rede receberá a condição inicial correspondente do vetor.
	* Se o número de elementos desses vectors for menor que o número de neurônios, a sequência de elementos do vector será repetida até que todos neurônios tenham tido um valor atribuído a eles.
	* Se o número de elementos desses vectors for zero, a distribuição será aleatória.
	* Os vectors x0 e y0 não precisam ter o mesmo número de elementos.
+ x0 e y0 podem ser omitidos para uma distribuição aleatória.


#### - Se a rede não for global, passar o vetor de Adjacência da Rede:
```c++
rede.readAdjVet( |AdjVet filePath| );
```

#### - (PLASTICIDADE) Se a rede for com plasticidade, defina os pesos de sinapse iniciais
```c++
rede.readSynWeights( |synWeights filePath| );
```


#### - (OPCIONAL) Para adicionar um sinal externo à rede, siga os passos.
1. Defina uma função void com o efeito do sinal em um dos neurônios afetados.
   Ela deve receber os seguintes parâmetros:
	+ double& x - referência do valor de x de um neurônio afetado no tempo n
	+ double& y - referência do valor de y de um neurônio afetado no tempo n
	+ int n 	 - tempo atual
	+ void* p  - ponteiro para um struct definido pelo usuário que contém os parâmetros do sinal. Deve ser convertido para o tipo apropriado antes de ser usado.

2.	Crie um objeto da classe SinalExterno e passe as características do sinal
	```c++
	SinalExterno sinal;
	sinal.setSinal(|Ponteiro da função definida acima|, [Ponteiro para o objeto que contém os parâmetros], [n de inicio do sinal]);
	sinal.setAfetados(|Num. de neurons na rede|, |Fração de afetados|, [Neurônios Afetados]);
	```
	+ "n de início do sinal" é o tempo a partir do qual o sinal será aplicado, acumulando os tempos de transiente e oficial.
	Zero inicia o sinal no primeiro tempo do transiente. Se o valor passado for igual ao transiente, os sinal é aplicado imediatamente após o fim do transiente.
	O valor default de "n de início do sinal" é zero.
	+ "neurônios afetados" é um vetor com os índices dos neurônios afetados. Pode ser omitido, ou passado como NULL para uma seleção aleatória.

3. Passe o endereço do objeto SinalExterno para a rede antes de chamar `rede.calcula()`.
	```c++
	rede.setSinal(&sinal);
	```

	Veja exemplo de código na seção correspondente para maior clareza.



#### - Defina as variáveis que serão calculadas e guardadas em memória no final da integração:
```c++
rede.setCalc( Grandeza::|Enum da Variável| );
```

Pode-se definir quantas variáveis forem necessárias e quaisquer umas.
Todas devem conter o prefixo 'Grandeza::'
Lista de Enums de Variáveis:

+ all 				- Todas
+ x 					- Variável x para cada tempo para cada neurônio
+ y 					- Variável y para cada tempo para cada neurônio
+ fase				- Fase de cada neurônio para todos os tempos
+ ordem 				- Parâmetro de Ordem da Rede para todos os tempos
+ campoMed			- Campo Médio da Rede para todos os tempos
+ varCampoMed			- Variância do Campo Médio
+ burstStart			- Tempos nos quais houve início de burst, para todos neurônios
+ onlyOrdemMed		- Calcula o Parâmetro de Ordem Médio sem guardar o Parâmetro de Ordem para todos os tempos
+ onlyAvgCampoMed		- Calcula a média temporal do Campo Médio sem guardar o Campo Médio para todos os tempos


#### - Mande calcular:
```c++
rede.calcula();
```



## Funcionalidades Adicionais


#### - Writters: Escrevem um arquivo .dat com a variável desejada.
+ Variáveis da Rede.
Exemplo: `rede.escreveOrdem( filePath, [header]);`

+ Variáveis do Neurônio.
Exemplo: `rede.escreveX( |filePath|, |Tamanho da Amostra|, [Neurônios escolhidos], [header] );`

Argumentos opicionais:
	Neurônios escolhidos: Vetor com os índices dos neuronios cuja variável será escrita. Default=NULL (Escolhe os primeiros)
	Header: Pode ser 0(No header), 1(Apenas header), 2(Header e Detalhes da Rede). Default=1


#### - Liberamento de Memória: Libera vetores alocados que não serão mais usados

Exemplos
```c++
rede.freeFase();
rede.destroy();	//Libera todos vetores
```

Não é obrigatória a chamada dessas funções para liberar a memória


#### - CleanUp: Se definido true, libera automaticamente os vetores auxiliares alocados
```c++
rede.setCleanUp( [bool] );
```

Default = false

[//]: <> (Exporta Rede: Salva os dados e vetores de uma rede em um arquivo .dat)
[//]: <> (	rede.exportRede( |filePath| );)

[//]: <> (Importa Rede: Lê os dados e vetores de uma rede salvos um arquivo .dat e os põe em um objeto rede)
[//]: <> (	rede.importRede( |filePath| );)

[//]: <> (	Ler os dados de um arquivo é aparentemente mais lento do que simplesmente integrar a rede.)
[//]: <> (	Esse feature é meio inútil =/)




# EXEMPLO DE CÓDIGOS

## TODO: Outros exemplos

#### - Exemplo de código com sinal externo:
```c++
//#includes

#define N 1000
#define SINAL_START_STEP 0			//Sinal começa a ser aplicado imediatamente no começo do transiente
#define FRACAO_AFETADOS_SINAL 0.2	//20% dos neurônios serão afetados

//Struct com os parâmetros do sinal (deve ser altera conforme o sinal desejado)
struct StrParametros {
	double frequencia, amplitude;
};

//Função de um sinal cossenoidal (deve ser altera conforme o sinal desejado)
void funcSinal(float &x, float &y, int n, void *p) {
	StrParametros *param = (StrParametros*) p;

	x += param->amplitude * cos( param->frequencia * n);
}

int main() {
	//Cria um struct com os parâmetros
	StrParametos param;
	param.frequencia = 0.001;
	param.amplitude = 0.01;

	//Escolhe os 20% primeiros neuronios para ser afetados (pode ser altera conforme o sinal desejado)
	int neuronsAfetados[FRACAO_AFETADOS_SINAL * N];
	for (int i=0; i < FRACAO_AFETADOS_SINAL * N; i++)
		neuronsAfetados[i] = i;

	//Cria o objeto SinalExterno
	SinalExterno sinal;
	sinal.setSinal(funcSinal, &param, SINAL_START_STEP);
	sinal.setAfetados(N, FRACAO_AFETADOS_SINAL, neuronsAfetados);

	//Cria rede
	Rede rede(|...|);
	rede.init();
	rede.setCalc(|...|);

	//Passa sinal para rede
	rede.setSinal(&sinal);

	rede.calcula

	return 0;
}
```


************************************************************************

