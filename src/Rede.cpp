#include <include/Rede.h>


//----- Constructors ------------------------------------------

Rede::Rede(uint numNeurons, uint numStep, uint transiente, double eps, double beta, double sig) :
RedeBase(numNeurons, numStep, transiente, eps, beta, sig) {

	adjVetSize = 0;
	conectividadeMedia = 0.0;
	adjVet = NULL;
	m_contribVizinhos = NULL;
}


Rede::Rede(uint numNeurons, uint numStep, uint transiente, double eps, const std::string &adjVetPath) :
RedeBase(numNeurons, numStep, transiente, eps) {

	adjVetSize = 0;
	conectividadeMedia = 0;
	adjVet = NULL;
	m_contribVizinhos = NULL;

	readAdjVet(adjVetPath);
}


Rede::~Rede() {
	destroy();
}


//----- Setters ----------------------------------------------

void Rede::readAdjVet(const std::string &filePath) {
    uint conexoes, neuronId;
    std::string auxStr;
    std::ifstream is(filePath.c_str());

    if (is.fail()) {
    	Log::Error("Failed to open adjVet " + filePath + "!");
    	exit(1);
    }

	//Get numero de linhas (elementos)
    adjVetSize = 0;
    while (getline(is, auxStr)) {
    	adjVetSize++;
    }

    is.clear();
    is.seekg(0, is.beg);

	//Aloca vetor
    if (adjVet != NULL)
        delete[] adjVet;
    adjVet =  new uint[adjVetSize]();
        
	conectividadeMedia = 0.0;
	conexoes = 0;
	neuronId = 0;

    //Enche vetor e acha a conectividade média
    for (uint i=0; i < adjVetSize; i++) {
        is >> adjVet[i];

    	//Checa por erros no ordem dos elementos de adjVet
        if (adjVet[i] - neuronId * numNeurons < 0) {
        	Log::Error("Erro no valor de adjVet no elemento %d!", i);
        	Log::Error("Checar arquivo " + filePath);
        	exit(1);
        }

    	//Checa se a ligação atual é do neurônio atual ou de um seguinte
        if (adjVet[i] < (neuronId + 1) * numNeurons) {
        	conexoes ++;
        } else {
            //Se já acabou as ligações do neurônio atual, adiciona o nº de conexões à soma pra média e reseta
            conectividadeMedia += conexoes;
            conexoes = 1;
            
            //Acha de qual neurônio é a ligação atual
            while (adjVet[i] < (neuronId + 1) * numNeurons)
        	   neuronId++;
        }
    }
    
    //Finaliza cálculo da conectividade média
    conectividadeMedia /= numNeurons;

    is.close();
}

void Rede::setAdjVet(uint *_adjVet, int _adjVetSize, double _conectividadeMedia) {
    adjVetSize = _adjVetSize;
    conectividadeMedia = _conectividadeMedia;
    
    adjVet = new uint[_adjVetSize];
    for (int i=0; i < _adjVetSize; ++i)
        adjVet[i] = _adjVet[i];
}


//----- Calculators ------------------------------------------

void Rede::calcula() {
    //----- Checagem de Erro -----
    bool noneSet = true;

    for (uint i=0; i < CALC_LIST_SIZE; i++) {
        if (m_shouldCalc[i] == true)
            noneSet = false;
    }
    //Se nada foi settado para ser calculado, termina
    if (noneSet) {
        Log::Error("Nenhuma grandeza foi escolhida para ser calculada!");
        exit(1);
    }

  	//Se não foi passado um vetor de adjacência, encerra
    if (adjVet == NULL) {
        Log::Error("Adj Vet não alocado!");
        exit(1);
    }

    //Checa se init() foi chamado
    if (m_alpha == NULL) {
        Log::Error("rede.init() não foi chamado!");
        exit(1);
    }
    //----------------------------


    	//Aloca vetor auxiliar m_contribVizinhos
    if (m_contribVizinhos == NULL)
    	m_contribVizinhos = new double[numNeurons];
    
    
    if (m_shouldCalc[Grandeza::none]) {
        for (uint k=0; k < CALC_LIST_SIZE; k++)
            m_shouldCalc[k] = false;
        
        m_shouldCalc[Grandeza::none] = true;
        
        calcNone();
    }

        //Se x e fase serão salvos
    if (m_shouldCalc[Grandeza::x] && m_shouldCalc[Grandeza::fase]) {
            //Calcula xy
        calcXY();

            //Se não quiser o y, libera sua memória
        if (!m_shouldCalc[Grandeza::y])
            freeYMat();

            //Calcula fase
        calcFaseDadoX();


            //Calcula Burst Starts
        if (m_shouldCalc[Grandeza::burstStart])
            calcBurstStartDadoX();

            //Calcula Parâmetro de Ordem
        if (m_shouldCalc[Grandeza::ordem])
            calcOrdemDadoFase();

            //Calcula Campo Médio
        if (m_shouldCalc[Grandeza::campoMed])
            calcCampoMedDadoX();
    }


        //Se x será salvo e fase não
    if (m_shouldCalc[Grandeza::x] && !m_shouldCalc[Grandeza::fase]) {
            //Calcula xy
        calcXY();

            //Se não quiser o y, libera sua memória
        if (!m_shouldCalc[Grandeza::y])
            freeYMat();

            //Calcula Burst Starts
        if (m_shouldCalc[Grandeza::burstStart])
            calcBurstStartDadoX();

            //Calcula Parâmetro de Ordem
        if (m_shouldCalc[Grandeza::ordem])
            calcOrdemDadoX();

            //Calcula Campo Médio
        if (m_shouldCalc[Grandeza::campoMed])
            calcCampoMedDadoX();
    }


        //Se fase será salva e x não
    if (!m_shouldCalc[Grandeza::x] && m_shouldCalc[Grandeza::fase]) {
            //Calcula Fase e Campo Médio
        if (m_shouldCalc[Grandeza::campoMed]) {
            calcFase_CampoMed();
        } else {
            calcFase();
        }

            //Calcula Parâmetro de Ordem
        if (m_shouldCalc[Grandeza::ordem])
            calcOrdemDadoFase();
    }


        //Se nem x nem fase serão salvos
    if (!m_shouldCalc[Grandeza::x] && !m_shouldCalc[Grandeza::fase]) {
            //Calcula Parâmetro de Ordem e Campo Médio
        if (m_shouldCalc[Grandeza::ordem] && m_shouldCalc[Grandeza::campoMed]) {
            calcOrdem_CampoMed();
        }
        else {
            if (m_shouldCalc[Grandeza::ordem])
                calcOrdem();

            if (m_shouldCalc[Grandeza::campoMed]) {
                if (m_shouldCalc[Grandeza::burstStart])
                    calcBurstStart_CampoMed();
                else
                    calcCampoMed();
            }

            if (m_shouldCalc[Grandeza::burstStart] && !m_shouldCalc[Grandeza::ordem] && !m_shouldCalc[Grandeza::campoMed])
                calcBurstStart();
        }
    }

        //Cálculo da Variância do Campo Médio
    if (m_shouldCalc[Grandeza::varCampoMed])
        calcVarCampoMed();


        //Limpa verores auxiliares da memória, caso não sejam mais usados
    if (m_cleanUp)
        freeAuxiliares();
}

void Rede::calcContribVizinhos_eletrico() {
	uint i, j, k;
	uint adjVetPos = 0;

    
    for (i=0; i < numNeurons; i++) {
        m_contribVizinhos[i] = 0;

        //Computa contribuição ao neuronio i dos vizinhos dele
        for (k = adjVetPos; k < adjVetSize; k++) {
            //Acha o próximo vizinho de i
            j = adjVet[k] - i * numNeurons;
            
            //Checa se j é mesmo um vizinho de i ou se é de (i+1)
            if (j >= numNeurons)
                break;
            
            //Calcula a contribuição de j em i
            m_contribVizinhos[i] += m_x[j];
        }

        adjVetPos = k;
    }


    for (i=0; i < numNeurons; i++)
    	m_contribVizinhos[i] /= conectividadeMedia;
}

void Rede::calcContribVizinhos() {
    uint i, j, k;
    uint adjVetPos = 0;
    
    for (i=0; i < numNeurons; i++) {
        m_contribVizinhos[i] = 0;

        //Computa contribuição ao neuronio i dos vizinhos dele
        for (k = adjVetPos; k < adjVetSize; k++) {
            //Acha o próximo vizinho de i
            j = adjVet[k] - i * numNeurons;
            
            //Checa se j é mesmo um vizinho de i ou se é de (i+1)
            if (j >= numNeurons)
                break;
            
            //--- Calcula a contribuição de j em i ---
            //Se j está disparando, injeta corrente
            if (m_x[j] > limiarSinapse)
                m_contribVizinhos[i] += (potencialSinapse - m_x[i]);
        }

        adjVetPos = k;
    }
}


void Rede::evoluiStep(uint n) {
	uint i;
	double tmp;

	//Calcula a contribuição dos vizinhos na dinâmica de cada neurônio e salva em m_contribVizinhos
	calcContribVizinhos();

        //Aplica a equação de movimento para cada neurônio
    for (i=0; i < numNeurons; i++) {
		tmp = m_x[i];

        m_x[i] = m_alpha[i] / (1.0 + m_x[i] * m_x[i]) + m_y[i]  +  (eps / conectividadeMedia) * m_contribVizinhos[i];
        m_y[i] = m_y[i] - sig * tmp - beta;
    }

        //Aplica o sinal externo se houver
    if (m_hasExternalSignal  &&  (int)n > ((int)m_sinal->startStep - (int)transiente))
        m_sinal->aplica(m_x, m_y, n);
    
    //Chama função externa do usuário
    externalFunction(n);
}


//----- Printers ---------------------------------------------

void Rede::escreveAdjMat(const std::string &filePath, bool asPGM) {
    uint i, j, k, vizinhoId;
    uint adjVetPos, adjMatPos;
    std::ofstream os(filePath.c_str());

    if (os.fail()) {
        Log::Warning("Não foi possível escrever adjMat! Failed to open file " + filePath + "!");
        return;
    }


    if (asPGM)
        os << "P2\n" << numNeurons << ' ' << numNeurons << "\n1\n";


    adjVetPos = 0;

    for (i=0; i < numNeurons; i++) {
        adjMatPos = 0;

        for (k = adjVetPos; k < adjVetSize; k++) {
            vizinhoId = adjVet[k] - i * numNeurons;
                
            if (vizinhoId >= numNeurons)
                break;


            for (j=adjMatPos; j < vizinhoId; j++)
                os << "0 ";

            os << "1 ";

            adjMatPos = j+1;
        }

            //Preenche do último 1 da linha até o final com 0
        for (j=adjMatPos; j < numNeurons; j++)
            os << "0 ";
        os << std::endl;        

            //Atualiza indice atual do adjVet
        adjVetPos = k;
    }

    os.close();
}


//----- Memoria ----------------------------------------------

void Rede::freeAuxiliares() {
    //Free m_x
    if (m_x != NULL) {
        delete[] m_x;
        m_x = NULL;
    }

        //Free m_y
    if (m_y != NULL) {
        delete[] m_y;
        m_y = NULL;
    }

        //Free m_real
    if (m_real != NULL) {
        delete[] m_real;
        m_real = NULL;
    }

        //Free m_imag
    if (m_imag != NULL) {
        delete[] m_imag;
        m_imag = NULL;
    }

        //Free m_burstId
    if (m_burstId != NULL) {
        delete[] m_burstId;
        m_burstId = NULL;
    }

        //Free m_histX
    if (m_histX != NULL) {
        for (uint i=0; i < numNeurons; i++) {
            if (m_histX[i] != NULL)
                delete[] m_histX[i];
        }

        delete[] m_histX;
        m_histX = NULL;
    }

    	//Free m_contribVizinhos
    if (m_contribVizinhos != NULL) {
    	delete[] m_contribVizinhos;
    	m_contribVizinhos = NULL;
    }
}

void Rede::destroy() {
    freeFase();
    freeOrdem();
    freeCampoMed();
    freeXYMat();
    freeBurstStart();
    freeAuxiliares();

        //Free Alpha
    if (m_alpha != NULL) {
        delete[] m_alpha;
        m_alpha = NULL;
    }

    if (adjVet != NULL) {
    	delete[] adjVet;
    	adjVet = NULL;
    }


    clearCalcList();
}

