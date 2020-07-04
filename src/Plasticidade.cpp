#include <include/Plasticidade.h>


//----- Constructors ------------------------------------------

RedePlasticidade::RedePlasticidade(uint numNeurons, uint numStep, uint transiente, double eps, enum Modelo::Modelo modelId) :
Rede(numNeurons, numStep, transiente, eps) {
    double _amplPotenc, _razaoDepPot, _tempoSaturacao;
    
    //----- Checagem de Erro -----
    if (modelId == Modelo::none) {
        Log::Error("Você deve passar o número identificador do modelo no constructor da RedePlasticidade!");
        Log::Info("As opções são: \n\t'Modelo::BTDP'   -   Modelo 1: Considera apenas os tempos de início de burst.\n\t'Modelo::Spike'   -   Modelo 2: Considera os spikes no meio dos bursts.");
        exit(1);
    }
    //----------------------------
    
    //Valores iniciais 
	limiarSinapse = 0.0;
    m_angCoefBTDP = 0.0;
    m_modelo = modelId;
    attXHistOnPlasticity = true;
	
    m_lastBurstStart = NULL;
	m_lastSpikesTime = NULL;
	synWeight.clear();
	m_neuronsBurstingNow.clear();
	
	//Valores padrão dos parâmetros
    if (modelId == Modelo::BTDP) {
        _amplPotenc  = 0.008;    //Maior do que o usado por Gjorgjieva (0.0005). Ajustado para não demorar tanto
        _razaoDepPot = 0.4;      //Verificado experimentalmente para retina
        _tempoSaturacao = 58;    //Calculado para ter razões próxima das de Butts de IBI e BurstLength
        
        setPlasticityParams(_amplPotenc, -_razaoDepPot*_amplPotenc, _tempoSaturacao);
    }
    
    if (modelId == Modelo::Spike)
        setPlasticityParams(0.0001, 0.00008, 20, 20);       
}

RedePlasticidade::~RedePlasticidade() {
	destroy();
}

void RedePlasticidade::destroy() {
	if (m_lastBurstStart != NULL) {
		delete[] m_lastBurstStart;
		m_lastBurstStart = NULL;
	}
    
    if (m_lastSpikesTime != NULL) {
        delete[] m_lastSpikesTime;
        m_lastSpikesTime = NULL;
    }
}


//-------------------- Setters --------------------
void RedePlasticidade::readSynWeights(const std::string &weightFilePath) {
    uint numElem;
    std::string auxStr;
    std::ifstream is(weightFilePath.c_str());

    if (is.fail()) {
        Log::Error("Failed to open adjVet " + weightFilePath + "!");
        exit(1);
    }

    //Get numero de linhas (elementos)
    numElem = 0;
    while (getline(is, auxStr)) {
        numElem++;
    }

    is.clear();
    is.seekg(0, is.beg);

    //----- Checagem erros -----
    if (adjVet == NULL) {
        Log::Error("O vetor de adjacência deve ser settado antes das forças de conexao!");
        is.close();
        exit(1);
    }
    
    if (numElem != adjVetSize) {
        Log::Error("O vector de força de conexões deve o mesmo número de elemtos que o vetor de adjacência!");
        is.close();
        exit(1);
    }
    //--------------------------
        
    
    synWeight.resize(adjVetSize);

    //Enche vetor e acha a conectividade média
    for (uint i=0; i < adjVetSize; i++)
        is >> synWeight[i];
    
    is.close();
    
    
    //Checagem de Erro
    for (uint i=0; i < adjVetSize; i++) {
        if (synWeight[i] < SYNWEIGHT_MIN  ||  synWeight[i] > SYNWEIGHT_MAX) {
            Log::Error("A força de conexão deve ser um valor entre %f e %f, porém o elemento %d passado tem valor %f!", SYNWEIGHT_MIN, SYNWEIGHT_MAX, i, synWeight[i]);
            exit(1);
        }
    }
}

void RedePlasticidade::setSynWeights(std::vector<double> _synWeight) {
	uint passedSize = _synWeight.size();

	//----- Checagem de Erro -----
	if (adjVet == NULL) {
		Log::Error("O vetor de adjacência deve ser settado antes das forças de conexao!");
		exit(1);
	}
	
	if (passedSize < 1) {
		Log::Error("O vector de força de conexões deve ter ao menos 1 elemento!");
		exit(1);
	}
	
	if (passedSize > adjVetSize) {
		Log::Warning("Foram passadas mais forças de conexões do que ha conexões! Usando apenas as primeiras.");
	}
	//----------------------------	
		
	//Se um valor foi passado
	if (passedSize == 1) {
		synWeight.resize(adjVetSize, _synWeight[0]);
		return;
	}
	
	//Se foram passados menos valores do que o necessário, repete os valores passados
	if (passedSize < adjVetSize) {
		synWeight.resize(adjVetSize);
		
		for (uint i=0; i < adjVetSize; i++)
			synWeight[i] = _synWeight[ i % passedSize ];
		
		return;
	}
	
	//Se foi passado um número igual ou maior que o necessário, preenche normalmente
	for (uint i=0; i < adjVetSize; i++) {
		synWeight.resize(adjVetSize);
	
		synWeight[i] = _synWeight[i];
	}
}

void RedePlasticidade::setPlasticityParams(double _amplPotenc, double _amplDepres, double _tauPotenc, double _tauDepres) {	
    //Checagem de Erro
    if (m_modelo != Modelo::Spike) {
        Log::Error("setPlasticityParams() só deve ser chamado com 4 argumentos se o modelo for Spike!");
        exit(1);
    }
    
    amplPotenc = _amplPotenc; 
	amplDepres = _amplDepres; 
	tauPotenc = _tauPotenc; 
	tauDepres = _tauDepres;
}

void RedePlasticidade::setPlasticityParams(double _amplPotenc, double _amplDepres, uint _tempoSaturacao) {
    //----- Checagem de Erro -----
    if (m_modelo != Modelo::BTDP) {
        Log::Error("setPlasticityParams() só deve ser chamado com 3 argumentos se o modelo for BTDP!");
        exit(1);
    }
    
    if (_amplDepres > 0.0)
        Log::Warning("Espera-se uma amplitude de depressão deve ser menor que zero, mas recebeu-se uma maior que zero! Continuando programa...");
    //----------------------------
    
    amplPotenc = _amplPotenc;
    amplDepres = _amplDepres;
    tempoSaturacao = _tempoSaturacao;
    
    m_amplDepInterna = amplDepres / 2.0;
    m_amplPotInterna = amplPotenc - m_amplDepInterna;
    m_angCoefBTDP = (m_amplPotInterna - m_amplDepInterna) / tempoSaturacao;
}


//-------------------- Calculators --------------------
void RedePlasticidade::calcula() {
    uint i;
    
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
        Log::Info("Chame o método readAdjVet().");
        exit(1);
    }
    
    //Se não foi passado um vetor de força de conexão, encerra
    if (synWeight.size() == 0) {
        Log::Error("Vetor de pesos de sinapses não alocado!");
        Log::Info("Chame o método setSynWeight().");
        exit(1);
    }

    //Checa se init() foi chamado
    if (m_alpha == NULL) {
        Log::Error("rede.init() não foi chamado!");
        exit(1);
    }
    //----------------------------


    //----- Alocação de Vetores Auxiliares -----
    if (m_contribVizinhos == NULL)
    	m_contribVizinhos = new double[numNeurons];
    
    if (m_modelo == Modelo::BTDP  &&  m_lastBurstStart == NULL) {
    	m_lastBurstStart = new int[numNeurons];
    	m_neuronsBurstingNow.reserve(numNeurons / 2);
        
        //Aloca m_histX, uma variável para o cálculo de início de burst definida em Base.h
        if (m_histX == NULL)
            m_histX = allocMat<float>(numNeurons, HIST_X_SIZE);
    }
    
    if (m_modelo == Modelo::Spike  &&  m_lastSpikesTime == NULL) {
        m_lastSpikesTime = new std::deque<uint>[numNeurons];
        for (i=0; i < numNeurons; i++)
            m_lastSpikesTime[i].resize(NUM_SPIKES_PLAS);
    }    
    //------------------------------------------
    

    //----- Inicialização de Variáveis -----
    if (m_modelo == Modelo::BTDP) {
        for (i=0; i < numNeurons; i++)
        	m_lastBurstStart[i] = 0;
        
        for (i=0; i < numNeurons; i++) {
            for (uint n=0; n < HIST_X_SIZE; n++)
                m_histX[i][n] = 10.0;
        }
    }
    
    if (m_modelo == Modelo::Spike) {
        for (i=0; i < numNeurons; i++) {
            for (uint k=0; k < NUM_SPIKES_PLAS; k++)
                m_lastSpikesTime[i][k] = 0;
        }
    }
    //--------------------------------------
    
        
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


void RedePlasticidade::evoluiStep(uint n) {
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
    
    
    //Atualiza força de conexões
    if (m_modelo == Modelo::BTDP)
        atualizaSynWeight_BTDPModel(n);
    if (m_modelo == Modelo::Spike)
        atualizaSynWeight_spikeModel(n);
    if (m_modelo == Modelo::Estatico) { }
        
    //Chama função externa do usuário
    externalFunction(n);
}

void RedePlasticidade::calcContribVizinhos() {
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
            	m_contribVizinhos[i] += synWeight[k] * (potencialSinapse - m_x[i]);
        }

		adjVetPos = k;
    }
}


void RedePlasticidade::tiraCITrans(int duracao) {
    //----- Checagem de Erro -----
    if (m_alpha == NULL) {
        Log::Error("init() deve ser chamada antes de tiraCITrans()!");
        exit(1);
    }
    
    //Se não foi passado um vetor de adjacência, encerra
    if (adjVet == NULL) {
        Log::Error("Adj Vet não alocado!");
        Log::Info("Chame o método readAdjVet() antes de tiraCITrans()");
        exit(1);
    }
    //----------------------------
    /*
    //Lê alpha atual e põe em um vector pra passar como argumento
    std::vector<double> alpha(numNeurons);
    for (uint i=0; i < numNeurons; ++i)
        alpha[i] = m_alpha[i];
    
    
    //Rede temporária que calcula o parâmetro de ordem antes da plasticidade agir
    Rede redeInicial(numNeurons, 1, duracao, 0.0, beta, sig);
    
    if (_eps > 0.0)
        redeInicial.setEps(_eps);
    
    redeInicial.setAdjVet(adjVet, adjVetSize, conectividadeMedia);
    redeInicial.init(alpha, m_initCond[0], m_initCond[1]);
    redeInicial.setCalc(Grandeza::none);
    redeInicial.calcula();
    
    m_initCond[0].resize(numNeurons);
    m_initCond[1].resize(numNeurons);
    
    //Atualiza condição inicial
    for (uint i=0; i < numNeurons; ++i) {
        m_initCond[0].at(i) = redeInicial.getAuxX()[i];
        m_initCond[1].at(i) = redeInicial.getAuxY()[i];
    }
    */
    
    Modelo::Modelo realModel = m_modelo;
    m_modelo = Modelo::Estatico;
    
    //----- Aloca Vetores -----
    if (m_x == NULL)
        m_x = new double[numNeurons];

    if (m_y == NULL)
        m_y = new double[numNeurons];
    
    if (m_contribVizinhos == NULL)
        m_contribVizinhos = new double[numNeurons];
    //-------------------------
    
    //--- Evolui rede ---
    m_isTrans = true;
    for (int n=0; n < duracao; n++)
        evoluiStep(n);
    m_isTrans = false;
    //-------------------
    
    m_initCond[0].resize(numNeurons);
    m_initCond[1].resize(numNeurons);
    
    //Atualiza condição inicial
    for (uint i=0; i < numNeurons; ++i) {
        m_initCond[0].at(i) = m_x[i];
        m_initCond[1].at(i) = m_y[i];
    }
    
    m_modelo = realModel;
}

//-------------------- Plasticidade -------------------
void RedePlasticidade::atualizaSynWeight_BTDPModel(uint n) {
	uint i, j, adjVetPos, idx_neuron;
    
    //Ajusta os tempos do último burst quando acaba o transiente
    if (!m_isTrans && n == 0) {
        for (i=0; i < numNeurons; i++)
            m_lastBurstStart[i] = m_lastBurstStart[i] - transiente;
    }
    
    //Acha neurônios que iniciaram burst no tempo atual
    m_neuronsBurstingNow.clear();
    for (i=0; i < numNeurons; i++) {
    	if (checaBurst(i))
    		m_neuronsBurstingNow.push_back(i);
    }
    
    /*Atualiza histX, a menos que ele já esteja sendo atualizado na função calc chamada; isto é
    quando o burstStart, a fase, ou a ordem estão sendo calculados, e x não. */
    if (attXHistOnPlasticity) {
        if (m_isTrans  ||  m_shouldCalc[Grandeza::x]  ||
        (!m_shouldCalc[Grandeza::burstStart] && !m_shouldCalc[Grandeza::fase] && !m_shouldCalc[Grandeza::ordem]) ) {
            for (i=0; i < numNeurons; i++)
                m_histX[i][ n % HIST_X_SIZE ] = m_x[i];
        }
    }
    
    
    //Não contabiliza no primeiro tempo, para evitar erros
    if (n == 0)
        return;

    
    //Para cada neurônio que iniciou burst nesse tempo
    for (i=0; i < m_neuronsBurstingNow.size(); i++) {
    	//Acha neurônio que está burstando
    	idx_neuron = m_neuronsBurstingNow[i];
    	
        //Atualiza o tempo de burst mais recente do neurônio atual
        m_lastBurstStart[idx_neuron] = n;
        
        
		//Acha a primeira sinapse que chega em idx_neuron (linha idx_neuron da adjMat)
		adjVetPos = 0;
		while (adjVet[adjVetPos] < idx_neuron*numNeurons)
			adjVetPos++;
		
		
		//----- Atualiza as sinapses que chegam em idx_neuron -----
        //Enquanto o adjVet atual for referente ao neurônio idx_neuron
		while (adjVet[adjVetPos] < (idx_neuron+1) * numNeurons  &&  adjVetPos < adjVetSize) {
			//Acha vizinho de idx_neuron
			j = adjVet[adjVetPos] - idx_neuron*numNeurons;
			
			//Atualiza força de sinapse
			synWeight[adjVetPos] += calcMudancaBTDP(n - m_lastBurstStart[j]);
			
			//Checa se a força da sinapse modificada foi abaixo do mínimo ou acima do máximo
			if (synWeight[adjVetPos] < SYNWEIGHT_MIN)
				synWeight[adjVetPos] = SYNWEIGHT_MIN;
            
            if (synWeight[adjVetPos] > SYNWEIGHT_MAX)
                synWeight[adjVetPos] = SYNWEIGHT_MAX;
			
			adjVetPos++;
		}
		
		//----- Atualiza as sinapses que saem de idx_neuron -----
        //Para cada sinapse
		for (adjVetPos=0; adjVetPos < adjVetSize; adjVetPos++) {
			//Checa se a sinapse atual sai de idx_neuron
			if (adjVet[adjVetPos] % numNeurons  ==  idx_neuron) {
				//Acha vizinho de idx_neuron
				j = adjVet[adjVetPos] / numNeurons;
				
				//Atualiza força de sinapse
				synWeight[adjVetPos] += calcMudancaBTDP(n - m_lastBurstStart[j]);
				
                //Checa se a força da sinapse modificada foi abaixo do mínimo ou acima do máximo
                if (synWeight[adjVetPos] < SYNWEIGHT_MIN)
                    synWeight[adjVetPos] = SYNWEIGHT_MIN;
                
                if (synWeight[adjVetPos] > SYNWEIGHT_MAX)
                    synWeight[adjVetPos] = SYNWEIGHT_MAX;
			}
		}
	}
}

void RedePlasticidade::atualizaSynWeight_spikeModel(uint n) { 
    uint i, j, k, adjVetPos;

    //Registra quais neurônios spikaram no tempo atual
    for (i=0; i < numNeurons; i++) {
        if (m_x[i] > LIMIAR_SPIKE) {
            m_lastSpikesTime[i].pop_front();  //Remove elemento mais antigo
            m_lastSpikesTime[i].push_back(n);
        }
    }
    
    //***** Atualiza as Conexões *****
        
    //Não contabiliza os primeiros tempos, para evitar erros
    if (n < 50)
        return;
    
    //Para cada neurônio
    for (i=0; i < numNeurons; i++) {
        //Se o neurônio não disparou no tempo atual, não faz nada
        if (m_lastSpikesTime[i].back() != n)
            continue;
        
        //Acha a primeira sinapse que sai de i
        adjVetPos = 0;
        while (adjVet[adjVetPos] < i*numNeurons)
            adjVetPos++;
        
        
        //Diminui a força das sinapses que saem de i
        while (adjVet[adjVetPos] < (i+1) * numNeurons  &&  adjVetPos < adjVetSize) {
            //Acha vizinho de i
            j = adjVet[adjVetPos] - i*numNeurons;
            
            //Para cada um dos NUM_SPIKES_PLAS últimos spikes de j, diminui o valor da sinapse
            for (k=0; k < NUM_SPIKES_PLAS; k++)
                synWeight[adjVetPos] -= amplDepres * exp( -(double)(n - m_lastSpikesTime[j].at(k)) / tauDepres );
            
            //Checa se a força da sinapse modificada foi abaixo do mínimo permitido
            if (synWeight[adjVetPos] < SYNWEIGHT_MIN)
                synWeight[adjVetPos] = SYNWEIGHT_MIN;
            
            adjVetPos++;
        }
        
        //Aumenta a força de sinapses que chegam em i
        for (adjVetPos=0; adjVetPos < adjVetSize; adjVetPos++) {
            //Checa se a sinapse atual chega em i
            if (adjVet[adjVetPos] % numNeurons  ==  i) {
                //Acha vizinho de i
                j = (adjVet[adjVetPos] - i) / numNeurons;
                
                //Para cada um dos NUM_SPIKES_PLAS últimos spikes de j, aumenta o valor da sinapse
                for (k=0; k < NUM_SPIKES_PLAS; k++)
                    synWeight[adjVetPos] += amplPotenc * exp( -(double)(n - m_lastSpikesTime[j].at(k)) / tauPotenc );
                
                //Checa se a força da sinapse modificada foi acima do máximo permitido
                if (synWeight[adjVetPos] > SYNWEIGHT_MAX)
                    synWeight[adjVetPos] = SYNWEIGHT_MAX;
            }
        }
    }   
}


double RedePlasticidade::calcMudancaBTDP(int burstStartLatency) {
    /* Aplica fórmula da curva BTDP */
    if (abs(burstStartLatency) > tempoSaturacao)
        return m_amplDepInterna;
    
    return m_amplPotInterna - fabs(m_angCoefBTDP * burstStartLatency);
}


//-------------------- Writters --------------------
void RedePlasticidade::escreveSynWeights(const std::string &filePath, int header, bool writeAdjVet) const {
    if (synWeight.size() == 0) {
        Log::Warning("Não é possível escrever Pesos de Sinapse. Vetor não alocado!");
        return;
    }

    std::ofstream outFile(filePath.c_str());

    if (outFile.fail()) {
        Log::Warning("Não foi possível abrir arquivo ", filePath.c_str());
        return;
    }

    //Imprime header
    if (header) {
        if (writeAdjVet)
            outFile << "adjVetIndex" << '\t' << "synWeight" << std::endl;
        else
            outFile << "synWeight" << std::endl;
    }

    //Imprime adjVet e pesos ao lado
    if (writeAdjVet) {
        for (uint k=0; k < adjVetSize; k++)
            outFile << adjVet[k] << '\t' << synWeight[k] << std::endl;;
    }
    //Imprime pesos
    else {
        for (uint k=0; k < adjVetSize; k++)
            outFile << synWeight[k] << std::endl;;
    }
    
    outFile.close();
}



