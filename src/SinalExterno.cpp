//#include <include/SinalExterno.h>
#include <include/SinalExterno.h>

//----- Constructors ----------------
SinalExterno::SinalExterno() {
	m_fracAfetados = 0.0;
	m_numAfetados = 0;
	startStep = 0;

	m_afetados = NULL;
	m_params = NULL;
	m_func = NULL;
}

SinalExterno::~SinalExterno() {
	if (m_afetados != NULL)
		delete[] m_afetados;
}


//----- Setters ---------------------
void SinalExterno::setSinal(SinalFuncPtr func, void *params, int startStep) {
	//Checagem de Erro
	if (startStep  < 0)
		startStep = 0;

	m_func = func;
	m_params = params;
	this->startStep = startStep;
}

void SinalExterno::setAfetados(uint numNeurons, double fracAfetados, int *chosenNeurons) {
	uint i, j;

	//Checagem de Erro
    if (fracAfetados > 1.0 || fracAfetados < 0.0) {
        Log::Error("Fração de neuronios afetados inválida: ", m_fracAfetados);
        exit(1);
    }


    //Salva informações na classe
	m_fracAfetados = fracAfetados;
	m_numAfetados = (uint)(fracAfetados * numNeurons);

	//(Re)Aloca vetor
    if (m_afetados != NULL)
    	delete[] m_afetados;
    m_afetados = new uint[m_numAfetados];



    if (chosenNeurons == NULL) {
            //Escolhe aleatoriamente os neuronios afetados pelo sinal
        m_afetados[0] = rand() % numNeurons;

        for (i=1; i < m_numAfetados; i++) {
            m_afetados[i] = rand() % numNeurons;

                //Checa se não foi escolhido um neuronio já afetado
            for (j=0; j < i; j++) {
                if (m_afetados[i] == m_afetados[j]) {
                    i--;
                    break;
                }
            }
        }
    }
    else {
        for (i=0; i < m_numAfetados; i++)
            m_afetados[i] = chosenNeurons[i];
    }
}


void SinalExterno::aplica(double *x, double *y, int n) {
	for (uint i=0; i < m_numAfetados; i++)
		m_func(x[ m_afetados[i] ], y [ m_afetados[i] ], n, m_params);
}