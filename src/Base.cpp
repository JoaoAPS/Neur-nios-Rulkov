#include <include/Base.h>


//-----Constructors------------------------------------------

RedeBase::RedeBase(uint numNeurons, uint numStep, uint transiente, double eps, double beta, double sig) :
numNeurons(numNeurons), numStep(numStep), transiente(transiente), eps(eps), beta(beta), sig(sig) {

    m_initCond.resize(2);
    ordemMed = 0;
    avgCampoMed = 0;
    varCampoMed = 0;
    alphaSeed = 0;
    ciSeed = 0;
    faseMin = 0;
    faseMax = 0;
    m_cleanUp = false;
    m_hasExternalSignal = false;
    m_inicioSinal = 0;

    m_alpha = NULL;
    ordem = NULL;
    fase = NULL;
    m_x = NULL;
    m_y = NULL;
    campoMed = NULL;
    xMat = NULL;
    yMat = NULL;
    m_sinal = NULL;
    m_burstId = NULL;
    m_histX = NULL;
    m_real = NULL;
    m_imag = NULL;
    burstStart = NULL;
    
    //Valores padrão
    setSynParams(LIMIAR_SPIKE, 1.0);

    clearCalcList();
}

RedeBase::~RedeBase() {
    destroy();
}

//-----Inits-------------------------------------------------

void RedeBase::init(std::vector<double> alpha, std::vector<double> x0, std::vector<double> y0) {
	//----- Checagem de Erro ------
	if (alpha.size() > numNeurons) {
        Log::Error("Tamanho do vector alpha inválido: ", alpha.size());
        exit(1);
	}

	if (x0.size() > numNeurons) {
        Log::Error("Tamanho do vector x0 inválido: ", x0.size());
        exit(1);
	}

	if (y0.size() > numNeurons) {
        Log::Error("Tamanho do vector y0 inválido: ", y0.size());
        exit(1);
	}
	//----------------------------

    //Aloca vetor alpha
    if (m_alpha == NULL)
        m_alpha = new double[numNeurons]();

    srand(alphaSeed);
    
    //Se o alpha passado tem tamanho zero, escolhe aleatório
    if (alpha.size() == 0) {
        for (uint i=0; i < numNeurons; i++)
            m_alpha[i] = 4.1 + (rand() % 3000000) / 10000000.0;
    }
    //Se tem numNeuron elementos, atribui os valores passado
    else if (alpha.size() == numNeurons) {
        for (uint i=0; i < numNeurons; i++)
            m_alpha[i] = alpha[i];
    }
    //Se apenas alguns valores de alpha foram passados
    else {
       for (uint i=0; i < numNeurons; i++)
            m_alpha[i] = alpha[ i % alpha.size() ];
    }
    
    
    //Condições iniciais
    m_initCond[0] = x0;
    m_initCond[1] = y0;
}

void RedeBase::init(float alpha, std::vector<double> x0, std::vector<double> y0) {
	init(std::vector<double>(1, alpha), x0, y0);
}

void RedeBase::init(std::vector<double> alpha, float x0, float y0) {
	init(alpha, std::vector<double>(1,x0), std::vector<double>(1,y0));
}

void RedeBase::init(float alpha, float x0, float y0) {
	init(std::vector<double>(1,alpha), std::vector<double>(1,x0), std::vector<double>(1,y0));
}

void RedeBase::clearCalcList() {
    for (uint i=0; i < CALC_LIST_SIZE; i++)
        m_shouldCalc[i] = false;
}


//----- Calculators ----------------------------------------------

void RedeBase::calcula() {
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

	//Checa se init() foi chamado
    if (m_alpha == NULL) {
        Log::Error("rede.init() não foi chamado!");
        exit(1);
    }
    //----------------------------
    
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


void RedeBase::calcNone() {
    uint n;

    //----- Aloca Vetores ----------------------------
    if (m_x == NULL)
        m_x = new double[numNeurons];

    if (m_y == NULL)
        m_y = new double[numNeurons];
    //-------------------------------------------------
    
       
    //Valores iniciais
    setInitialCondition();

    //Calcula transiente
    calcTransiente();

    //Evolui rede
    for (n=0; n < numStep; n++)
        evoluiStep(n);
}

void RedeBase::calcXY() {
    uint i, n;

    //----- Aloca Vetores -------------------------------
            //Vetor x auxiliar
    if (m_x == NULL)
        m_x = new double[numNeurons];

        //Vetor y auxiliar
    if (m_y == NULL)
        m_y = new double[numNeurons];


    if (xMat == NULL)
        xMat =  allocMat<float>(numNeurons, numStep);

    if (yMat == NULL)
        yMat =  allocMat<float>(numNeurons, numStep);
    //---------------------------------------------------

        //Condições iniciais
    setInitialCondition();


        //Calcula transiente
    calcTransiente();

 
     //Para cada step
    for (n=0; n < numStep; n++) {
        //Passa os resultados para as matrizes definitivas
        for (i=0; i < numNeurons; i++) {
            xMat[i][n] = m_x[i];
            yMat[i][n] = m_y[i];
        }
        
        //Evolui m_x e m_y em um 
        evoluiStep(n);
    }
}

void RedeBase::calcBurstStart() {
    uint i, n;

    //----- Aloca Vetores -----
    if (m_x == NULL)
        m_x = new double[numNeurons];

    if (m_y == NULL)
        m_y = new double[numNeurons];

    if (m_histX == NULL)
        m_histX = allocMat<float>(numNeurons, HIST_X_SIZE);

    if (burstStart == NULL)
        burstStart = new std::vector<uint>[numNeurons];
    //-------------------------


    //----- Valoroes inicias --------------------
    for (i=0; i < numNeurons; i++) {
        for (n=0; n < HIST_X_SIZE; n++)
            m_histX[i][n] = 10.0;

        //Reserva espaço mínimo estimado (300 é próximo do maior período de burst já visto)
        burstStart[i].clear();
        burstStart[i].reserve(numStep / 300);
    }

    faseMin = 0;
    faseMax = numStep - 1;
    
    setInitialCondition();
    //-------------------------------------------


        //Calcula transiente
    calcTransiente();


    //Para cada step
    for (n=0; n < numStep; n++) {
        //Evolui m_x e m_y em um step
        evoluiStep(n);
        

        //Para cada neurônio
        for (i=0; i < numNeurons; i++) {
            //Se o tempo atual é inicio de burst, adiciona à lista
            if  ( checaBurst(i) )
                burstStart[i].push_back(n);

            //Atualiza histX
            m_histX[i][ n % HIST_X_SIZE ] = m_x[i];
        }
    }//step
}

void RedeBase::calcFase() {
    uint i, n, k;

    //----- Aloca Vetores -----------------------
    allocForFase();

    if (fase == NULL)
        fase =  allocMat<float>(numNeurons, numStep);

    if (m_shouldCalc[Grandeza::burstStart] && burstStart == NULL)
        burstStart = new std::vector<uint>[numNeurons];
    //-------------------------------------------


    //----- Valoroes inicias --------------------
    for (i=0; i < numNeurons; i++) {
        m_burstId[i] = 0;

        for (n=0; n < HIST_X_SIZE; n++)
            m_histX[i][n] = 10.0;
    }

    if (m_shouldCalc[Grandeza::burstStart]) {
        //Reserva espaço mínimo estimado (300 é próximo do maior período de burst já visto)
        for (i=0; i < numNeurons; i++) {
            burstStart[i].clear();
            burstStart[i].reserve(numStep / 300);
        }
    }

    faseMin = 0;
    faseMax = numStep - 1;

    setInitialCondition();
    //-------------------------------------------


        //Calcula transiente
    calcTransiente();


        //Para cada step
    for (n=0; n < numStep; n++) {
            //Evolui m_x e m_y em um step
        evoluiStep(n);
        

            //Para cada neurônio
        for (i=0; i < numNeurons; i++) {
                //Se o tempo atual é inicio de burst
            if  ( checaBurst(i) ) {
                //Guarda o índice de burstStart, se necessário
                if (m_shouldCalc[Grandeza::burstStart])
                    burstStart[i].push_back(n);

                    //Checa se é o primeiro burst e enche de zeros se for
                if (m_burstId[i] == 0) {
                    for(k=0; k <= n; k++)
                        fase[i][k] = 0.0;

                        //Checa se é o primeiro burst mais atrasado
                    if (n > faseMin)
                        faseMin = n;
                }
                else {
                        //Se não for o primeiro, interpola a fase
                    for(k = m_burstId[i]; k <= n; k++) {
                        fase[i][k] =  fase[i][m_burstId[i]] + 2*PI * (k - m_burstId[i]) / (float)(n - m_burstId[i]);
                    }
                }

                    //Atualiza o indice do ultimo burst
                m_burstId[i] = n;
            }

                //Atualiza histX
            m_histX[i][ n % HIST_X_SIZE ] = m_x[i];
        }
    }//step


    for (i=0; i < numNeurons; i++) {
            //Acha primeiro último_burst
        if (m_burstId[i] < faseMax)
            faseMax = m_burstId[i];

            //Completa vetor fase
        for (n = m_burstId[i]; n < numStep - 1; n++) {
            fase[i][n+1] = fase[i][n];
        }
    }
}

void RedeBase::calcOrdem() {
    uint i, n, k;
    double fase, sum;

    //----- Aloca vetores -----
    allocForFase();
    allocForOrdem();

    if (m_shouldCalc[Grandeza::burstStart] && burstStart == NULL)
        burstStart = new std::vector<uint>[numNeurons];
    //-------------------------


    //----- Valores Iniciais ---------------------
    for (n=0 ; n < numStep; n++) {   
        m_real[n] = 0.0;
        m_imag[n] = 0.0;
    }

    for (i=0; i < numNeurons; i++) {
        m_burstId[i] = 0;

        for (n=0; n < HIST_X_SIZE; n++)
            m_histX[i][n] = 10.0;
    }

    if (m_shouldCalc[Grandeza::burstStart]) {
        //Reserva espaço mínimo estimado (300 é próximo do maior período de burst já visto)
        for (i=0; i < numNeurons; i++){
            burstStart[i].clear();
            burstStart[i].reserve(numStep / 300);
        }
    }

    faseMin = 0;
    faseMax = numStep - 1;

    setInitialCondition();
    //--------------------------------------------


        //Calcula transiente
    calcTransiente();


        //Para cada step
    for (n=0; n < numStep; n++) {
            //Evolui m_x e m_y em um step
        evoluiStep(n);
        
        
            //Para cada neurônio
        for (i=0; i < numNeurons; i++) {               
                //Se o tempo atual é inicio de burst
            if ( checaBurst(i) ) {
                //Guarda o índice de burstStart, se necessário
                if (m_shouldCalc[Grandeza::burstStart])
                    burstStart[i].push_back(n);

                    //Checa se é o primeiro burst
                if (m_burstId[i] == 0) {
                        //Checa se é o primeiro burst mais atrasado
                    if (n > faseMin)
                        faseMin = n;
                }
                else {
                    //Se não for o primeiro
                    for(k = m_burstId[i]; k < n; k++) {
                            //Interpola a fase para o tempo k
                        fase = 2*PI * (k - m_burstId[i]) / (float)(n - m_burstId[i]);
                        
                            //Adiciona às componenentes de ordem[k]
                        m_real[k] += cos(fase);
                        m_imag[k] += sin(fase);
                    }
                }

                    //Atualiza o indice do ultimo burst
                m_burstId[i] = n;
            }

                //Atualiza histX
            m_histX[i][ n % HIST_X_SIZE ] = m_x[i];
        }
    }//step


        //Acha primeiro último_burst
    for (i=0; i < numNeurons; i++) {
        if (m_burstId[i] < faseMax)
            faseMax = m_burstId[i];
    }


    sum = 0.0;

    if (m_shouldCalc[Grandeza::onlyOrdemMed]) {
        for (n = faseMin; n < faseMax; n++)
            sum += sqrt(pow(m_real[n], 2) + pow(m_imag[n], 2)) / numNeurons;
    }
    else {
            //Aloca ordem vet
        if (ordem == NULL)
            ordem = new double[numStep];

            //Antes de faseMin, ordem = 0
        for (n=0; n < faseMin; n++)
            ordem[n] = 0.0;

            //Depois de faseMax, ordem = 0
        for (n = faseMax; n < numStep; n++)
            ordem[n] = 0.0;

            //No intervalo válido, calcula a ordem
        for (n = faseMin; n < faseMax; n++) {
            ordem[n] = sqrt(pow(m_real[n], 2) + pow(m_imag[n], 2)) / numNeurons;

            sum += ordem[n];
        }
    }

        //Calcula Ordem Médio
    ordemMed = sum / (faseMax - faseMin);
}

void RedeBase::calcCampoMed() {
    uint i, n;
    double media, campoMedSum;


    //----- Aloca Vetores ----------------------------
    if (!m_shouldCalc[Grandeza::onlyAvgCampoMed] && campoMed == NULL)
        campoMed = new double[numStep];

    if (m_x == NULL)
        m_x = new double[numNeurons];

    if (m_y == NULL)
        m_y = new double[numNeurons];
    //-------------------------------------------------

    //Valores iniciais
    campoMedSum = 0.0;
    setInitialCondition();


        //Calcula transiente
    calcTransiente();


        //Para cada step
    for (n=0; n < numStep; n++) {
            //Evolui m_x e m_y em um step
        evoluiStep(n);

            //Calcula Campo Médio
        media = 0.0;
        for (i=0; i < numNeurons; i++)
            media += m_x[i];

        media /= numNeurons;

        if (!m_shouldCalc[Grandeza::onlyAvgCampoMed])
            campoMed[n] = media;

        campoMedSum += media;
    }

        //Calcula média temporal do Campo Médio
    avgCampoMed = campoMedSum / numStep;
}

void RedeBase::calcBurstStart_CampoMed() {
    uint i, n;
    double media, campoMedSum;

    //----- Aloca Vetores -----
    if (m_x == NULL)
        m_x = new double[numNeurons];

    if (m_y == NULL)
        m_y = new double[numNeurons];

    if (m_histX == NULL)
        m_histX = allocMat<float>(numNeurons, HIST_X_SIZE);

    if (burstStart == NULL)
        burstStart = new std::vector<uint>[numNeurons];

    if (!m_shouldCalc[Grandeza::onlyAvgCampoMed] && campoMed == NULL)
        campoMed = new double[numStep];
    //-------------------------


    //----- Valoroes inicias --------------------
    for (i=0; i < numNeurons; i++) {
        for (n=0; n < HIST_X_SIZE; n++)
            m_histX[i][n] = 10.0;

        //Reserva espaço mínimo estimado (300 é próximo do maior período de burst já visto)
        burstStart[i].clear();
        burstStart[i].reserve(numStep / 300);
    }

    faseMin = 0;
    faseMax = numStep - 1;

    campoMedSum = 0.0;
    
    setInitialCondition();
    //-------------------------------------------


        //Calcula transiente
    calcTransiente();


    //Para cada step
    for (n=0; n < numStep; n++) {
            //Evolui m_x e m_y em um step
        evoluiStep(n);

        //----- Calcula Campo Médio -----
        media = 0.0;
        for (i=0; i < numNeurons; i++)
            media += m_x[i];

        media /= numNeurons;

        if (!m_shouldCalc[Grandeza::onlyAvgCampoMed])
            campoMed[n] = media;

        campoMedSum += media;
        //-------------------------------
        

        //----- Calcula Burst Starts -----
        for (i=0; i < numNeurons; i++) {
                //Se o tempo atual é inicio de burst, adiciona à lista
            if  ( checaBurst(i) )
                burstStart[i].push_back(n);

                //Atualiza histX
            m_histX[i][ n % HIST_X_SIZE ] = m_x[i];
        }
        //--------------------------------
    }//step

        //Calcula média temporal do Campo Médio
    avgCampoMed = campoMedSum / numStep;
}

void RedeBase::calcFase_CampoMed() {
    uint i, n, k;
    double media, campoMedSum;

    //----- Aloca Vetores -----------------------
    allocForFase();

    if (fase == NULL)
        fase =  allocMat<float>(numNeurons, numStep);

    if (m_shouldCalc[Grandeza::burstStart] && burstStart == NULL)
        burstStart = new std::vector<uint>[numNeurons];
    //-------------------------------------------


    //----- Valores inicias --------------------
    for (i=0; i < numNeurons; i++) {
        m_burstId[i] = 0;

        for (n=0; n < HIST_X_SIZE; n++)
            m_histX[i][n] = 10.0;
    }

    faseMin = 0;
    faseMax = numStep - 1;

    if (m_shouldCalc[Grandeza::burstStart]) {
        //Reserva espaço mínimo estimado (300 é próximo do maior período de burst já visto)
        for (i=0; i < numNeurons; i++){
            burstStart[i].clear();
            burstStart[i].reserve(numStep / 300);
        }
    }

    setInitialCondition();
    //-------------------------------------------


        //Calcula transiente
    calcTransiente();


    campoMedSum = 0.0;

        //Para cada step
    for (n=0; n < numStep; n++) {
            //Evolui m_x e m_y em um step
        evoluiStep(n);


        //----- Calcula Campo Médio do tempo atual ---------
        media = 0.0;
        for (i=0; i < numNeurons; i++)
            media += m_x[i];

        media /= numNeurons;

        if (!m_shouldCalc[Grandeza::onlyAvgCampoMed])
            campoMed[n] = media;

        campoMedSum += media;
        //--------------------------------------------------


            //Para cada neurônio
        for (i=0; i < numNeurons; i++) {
                //Se o tempo atual é inicio de burst
            if  ( checaBurst(i) ) {
                //Guarda o índice de burstStart, se necessário
                if (m_shouldCalc[Grandeza::burstStart])
                    burstStart[i].push_back(n);

                    //Checa se é o primeiro burst e enche de zeros se for
                if (m_burstId[i] == 0) {
                    for(k=0; k <= n; k++)
                        fase[i][k] = 0.0;

                        //Checa se é o primeiro burst mais atrasado
                    if (n > faseMin)
                        faseMin = n;
                }
                else {
                        //Se não for o primeiro, interpola a fase
                    for(k = m_burstId[i]; k <= n; k++) {
                        fase[i][k] =  fase[i][m_burstId[i]] + 2*PI * (k - m_burstId[i]) / (double)(n - m_burstId[i]);
                    }
                }

                    //Atualiza o indice do ultimo burst
                m_burstId[i] = n;
            }

                //Atualiza histX
            m_histX[i][ n % HIST_X_SIZE ] = m_x[i];
        }
    }//step


        //Calcula média temporal do Campo Médio
    avgCampoMed = campoMedSum / numStep;


    for (i=0; i < numNeurons; i++) {
            //Acha primeiro último_burst
        if (m_burstId[i] < faseMax)
            faseMax = m_burstId[i];

            //Completa vetor fase
        for (n = m_burstId[i]; n < numStep - 1; n++) {
            fase[i][n+1] = fase[i][n];
        }
    }
}

void RedeBase::calcOrdem_CampoMed() {
    uint i, n, k;
    double fase, sum, media, campoMedSum;

    //----- Aloca vetores -----
    allocForFase();
    allocForOrdem();

    if (m_shouldCalc[Grandeza::burstStart] && burstStart == NULL)
        burstStart = new std::vector<uint>[numNeurons];
    
    if (!m_shouldCalc[Grandeza::onlyAvgCampoMed] && campoMed == NULL)
        campoMed = new double[numStep];
    //--------------------------

    //----- Valores Iniciais ---------------------
    for (n=0 ; n < numStep; n++) {   
        m_real[n] = 0.0;
        m_imag[n] = 0.0;
    }

    for (i=0; i < numNeurons; i++) {
        m_burstId[i] = 0;

        for (n=0; n < HIST_X_SIZE; n++)
            m_histX[i][n] = 10.0;
    }
    
    if (m_shouldCalc[Grandeza::burstStart]) {
        //Reserva espaço mínimo estimado (300 é próximo do maior período de burst já visto)
        for (i=0; i < numNeurons; i++){
            burstStart[i].clear();
            burstStart[i].reserve(numStep / 300);
        }
    }

    faseMin = 0;
    faseMax = numStep - 1;

    setInitialCondition();
    //--------------------------------------------


        //Calcula transiente
    calcTransiente();


    campoMedSum = 0.0;

        //Para cada step
    for (n=0; n < numStep; n++) {
            //Evolui m_x e m_y em um step
        evoluiStep(n);


        //----- Calcula Campo Médio do tempo atual ---------
        media = 0.0;
        for (i=0; i < numNeurons; i++)
            media += m_x[i];

        media /= numNeurons;

        if (!m_shouldCalc[Grandeza::onlyAvgCampoMed])
            campoMed[n] = media;

        campoMedSum += media;
        //--------------------------------------------------


            //Para cada neurônio
        for (i=0; i < numNeurons; i++) {
                //Se o tempo atual é inicio de burst
            if ( checaBurst(i) ) {
                //Guarda o índice de burstStart, se necessário
                if (m_shouldCalc[Grandeza::burstStart])
                    burstStart[i].push_back(n);

                    //Checa se é o primeiro burst
                if (m_burstId[i] == 0) {
                        //Checa se é o primeiro burst mais atrasado
                    if (n > faseMin)
                        faseMin = n;
                }
                else {
                    //Se não for o primeiro
                    for(k = m_burstId[i]; k < n; k++) {
                            //Interpola a fase para o tempo k
                        fase = 2*PI * (k - m_burstId[i]) / (double)(n - m_burstId[i]);
                        
                            //Adiciona às componenentes de ordem[k]
                        m_real[k] += cos(fase);
                        m_imag[k] += sin(fase);
                    }
                }

                    //Atualiza o indice do ultimo burst
                m_burstId[i] = n;
            }

                //Atualiza histX
            m_histX[i][ n % HIST_X_SIZE ] = m_x[i];
        }
    }//step


        //Calcula média temporal do Campo Médio
    avgCampoMed = campoMedSum / numStep;

        //Acha primeiro último_burst
    for (i=0; i < numNeurons; i++) {
        if (m_burstId[i] < faseMax)
            faseMax = m_burstId[i];
    }


    sum = 0.0;

    if (m_shouldCalc[Grandeza::onlyOrdemMed]) {
        for (n = faseMin; n < faseMax; n++)
            sum += sqrt(pow(m_real[n], 2) + pow(m_imag[n], 2)) / numNeurons;
    }
    else {
            //Aloca ordem vet
        if (ordem == NULL)
            ordem = new double[numStep];

            //Antes de faseMin, ordem = 0
        for (n=0; n < faseMin; n++)
            ordem[n] = 0.0;

            //Depois de faseMax, ordem = 0
        for (n = faseMax; n < numStep; n++)
            ordem[n] = 0.0;

            //No intervalo válido, calcula a ordem
        for (n = faseMin; n < faseMax; n++) {
            ordem[n] = sqrt(pow(m_real[n], 2) + pow(m_imag[n], 2)) / numNeurons;

            sum += ordem[n];
        }
    }

        //Calcula Ordem Médio
    ordemMed = sum / (faseMax - faseMin);
}

void RedeBase::calcBurstStartDadoX() {
    uint i, n, k;
    bool isBurst;

        //Checa se x foi dado
    if (xMat == NULL) {
        Log::Error("Não é possível calcular burstStart! xMat não alocada!");
        exit(1);
    }

    //----- Aloca Vetores -----
    if (burstStart == NULL)
        burstStart = new std::vector<uint>[numNeurons];
    //-------------------------


    //----- Valores Inicias -----
    for (i=0; i < numNeurons; i++) {
        //Reserva espaço mínimo estimado (300 é próximo do maior período de burst já visto)
        burstStart[i].clear();
        burstStart[i].reserve(numStep / 300);
    }
    //---------------------------


        //Para cada neurônio
    for (i=0; i < numNeurons; i++) {
            //Para cada step a partir do menor n que pode ser burst
        for (n = HIST_X_SIZE; n < numStep; n++) {

                //Checa se o tempo atual é o inicio de um burst
            isBurst = false;
            if (xMat[i][n] > LIMIAR_SPIKE) {
                isBurst = true;

                for (k = n - HIST_X_SIZE; k < n; k++) {
                    if (xMat[i][k] > LIMIAR_SPIKE) {
                        isBurst = false;
                        break;
                    }
                }
            }

            if (isBurst)
                burstStart[i].push_back(n);
        }
    }
}

void RedeBase::calcFaseDadoX() {
    uint i, n, k, burstId;
    bool isBurst;

        //Checa se x foi dado
    if (xMat == NULL) {
        Log::Error("Não é possível calcular fase! xMat não alocada!");
        exit(1);
    }

    //----- Aloca Vetores -----
    if (fase == NULL)
    	fase = allocMat<float>(numNeurons, numStep);
    //-------------------------


    //----- Valores Inicias -----
    faseMin = 0;
    faseMax = numStep - 1;
    //---------------------------


        //Para cada neurônio
    for (i=0; i < numNeurons; i++) {
        burstId = 0;

            //Para cada step a partir do menor n que pode ser burst
        for (n = HIST_X_SIZE; n < numStep; n++) {

                //Checa se o tempo atual é o inicio de um burst
            isBurst = false;
            if (xMat[i][n] > LIMIAR_SPIKE) {
                isBurst = true;

                for (k = n - HIST_X_SIZE; k < n; k++) {
                    if (xMat[i][k] > LIMIAR_SPIKE) {
                        isBurst = false;
                        break;
                    }
                }
            }

                //Se é inicio de burst
            if  (isBurst) {
                    //Checa se é o primeiro burst e enche de zeros se for
                if (burstId == 0) {
                    for(k=0; k <= n; k++)
                        fase[i][k] = 0.0;

                        //Checa se é o primeiro burst mais atrasado
                    if (n > faseMin)
                        faseMin = n;
                }
                else {
                        //Se não for o primeiro, interpola a fase
                    for(k = burstId; k <= n; k++) {
                        fase[i][k] =  fase[i][burstId] + 2*PI * (k - burstId) / (double)(n - burstId);
                    }
                }

                    //Atualiza o indice do ultimo burst
                burstId = n;

                    //Não haverá um novo burst nesse intervalo, então não calcula ele
                n += HIST_X_SIZE;
            }
        }//step

            //Checa se é o primeiro último_burst
        if (burstId < faseMax)
            faseMax = burstId;
            
            //Completa vetor fase
        for (n = burstId; n < numStep - 1; n++) {
            fase[i][n+1] = fase[i][n];
        }
    }//neuron
}

void RedeBase::calcOrdemDadoFase() {
    uint i, n;
    double real, imag;
    double sum = 0.0;

    //Checa se fase foi dado
    if (fase == NULL) {
        Log::Error("Não é possível calcular ordem! fase não alocado!");
        exit(1);
    }

    if (m_shouldCalc[Grandeza::onlyOrdemMed]) {
        for (n = faseMin; n < faseMax; n++) {
            real = 0.0;
            imag = 0.0;

                //Calcula parte real e imaginária de e^(i phi)
            for (i=0; i < numNeurons; i++) {
                real += cos(fase[i][n]);
                imag += sin(fase[i][n]);
            }
            
            sum += sqrt(pow(real, 2) + pow(imag, 2)) / numNeurons;
        }
    }

    else {
        if (ordem == NULL)
            ordem = new double[numStep]();

        for (n = faseMin; n < faseMax; n++) {
            real = 0.0;
            imag = 0.0;

                //Calcula parte real e imaginária de e^(i phi)
            for (i=0; i < numNeurons; i++) {
                real += cos(fase[i][n]);
                imag += sin(fase[i][n]);
            }

                //Calcula o módulo da soma de todos fasores no tempo j, e divide por numNeurons para deixar entre 0 e 1
            ordem[n] = sqrt(pow(real, 2) + pow(imag, 2)) / numNeurons;

            sum += ordem[n];
        }
    }

        //Calcula parâmetro de ordem médio no tempo
    ordemMed = sum / (faseMax - faseMin);
}

void RedeBase::calcOrdemDadoX() {
    uint i, n, k, burstId;
    double fase, sum;
    bool isBurst;

    //Checa se x foi dado
    if (xMat == NULL) {
        Log::Error("Não é possível calcular ordem! xMat não alocada!");
        exit(1);
    }

    //----- Aloca Vetores -----
    allocForOrdem();
    //-------------------------


    //----- Valores Iniciais ---------------------
    for (n=0 ; n < numStep; n++) {   
        m_real[n] = 0.0;
        m_imag[n] = 0.0;
    }

    faseMin = 0;
    faseMax = numStep - 1;
    //--------------------------------------------

        //Para cada neurônio
    for (i=0; i < numNeurons; i++) {
        burstId = 0;

            //Para cada step a partir do menor n que pode ser burst
        for (n = HIST_X_SIZE; n < numStep; n++) {

                //Checa se o tempo atual é o inicio de um burst
            isBurst = false;
            if (xMat[i][n] > LIMIAR_SPIKE) {
                isBurst = true;

                for (k = n - HIST_X_SIZE; k < n; k++) {
                    if (xMat[i][k] > LIMIAR_SPIKE) {
                        isBurst = false;
                        break;
                    }
                }
            }

                //Se é inicio de burst
            if (isBurst) {
                    //Checa se é o primeiro burst e enche de zeros se for
                if (burstId == 0) {
                        //Checa se é o primeiro burst mais atrasado
                    if (n > faseMin)
                        faseMin = n;
                }
                else {
                        //Se não for o primeiro, interpola a fase
                    for(k = burstId; k < n; k++) {
                        fase = 2*PI * (k - burstId) / (double)(n - burstId);

                            //Adiciona às componenentes de ordem[k]
                        m_real[k] += cos(fase);
                        m_imag[k] += sin(fase);
                    }
                }

                    //Atualiza o indice do ultimo burst
                burstId = n;

                    //Não haverá um novo burst nesse intervalo, então não calcula ele
                n += HIST_X_SIZE;
            }
        }//step

            //Checa se é o primeiro último_burst
        if (burstId < faseMax)
            faseMax = burstId;
    }


    sum = 0.0;

    if (m_shouldCalc[Grandeza::onlyOrdemMed]) {
        for (n = faseMin; n < faseMax; n++)
            sum += sqrt(pow(m_real[n], 2) + pow(m_imag[n], 2)) / numNeurons;
    }
    else {
            //Aloca ordem vet
        if (ordem == NULL)
            ordem = new double[numStep];

            //Antes de faseMin, ordem = 0
        for (n=0; n < faseMin; n++)
            ordem[n] = 0.0;

            //Depois de faseMax, ordem = 0
        for (n = faseMax; n < numStep; n++)
            ordem[n] = 0.0;

            //No intervalo válido, calcula a ordem
        for (n = faseMin; n < faseMax; n++) {
            ordem[n] = sqrt(pow(m_real[n], 2) + pow(m_imag[n], 2)) / numNeurons;

            sum += ordem[n];
        }
    }

        //Calcula Ordem Médio
    ordemMed = sum / (faseMax - faseMin);
}

void RedeBase::calcCampoMedDadoX() {
    uint i, n;
    avgCampoMed = 0.0;

            //Checa se x foi dado
    if (xMat == NULL) {
        Log::Error("Não é possível calcular campoMed! xMat não alocada!");
        exit(1);
    }


    if (m_shouldCalc[Grandeza::onlyAvgCampoMed]) {
        for (i=0; i < numNeurons; i++) {
            for (n=0; n < numStep; n++) {
                avgCampoMed += xMat[i][n];
            }
        }
        
        avgCampoMed /= (numNeurons * numStep);
    }

    else {
            //Aloca caso não já
        if (campoMed == NULL)
            campoMed = new double[numStep];

            //Zera vetor
        for (n=0; n < numStep; n++)
            campoMed[n] = 0;

            //Calcula 
        for (n=0; n < numStep; n++) {
            for (i=0; i < numNeurons; i++) {
                campoMed[n] += xMat[i][n];
            }

            campoMed[n] /= numNeurons;

            avgCampoMed += campoMed[n];
        }

        avgCampoMed /= numStep;
    }
}

void RedeBase::calcVarCampoMed() {
    if (campoMed == NULL) {
        Log::Error("Não é possível calcular variancia do campoMed! Vetor não alocado!");
        return;
    }

    double sum = 0.0;

    for (uint n=0; n < numStep; n++) {
        sum += (campoMed[n] - avgCampoMed) * (campoMed[n] - avgCampoMed);
    }

    varCampoMed = sum / numStep;
}


void RedeBase::calcTransiente() {
    m_isTrans = true;
    
	for (uint n=0; n < transiente; n++)
		evoluiStep(n);

	if (m_hasExternalSignal)
		m_inicioSinal = m_sinal->startStep - transiente;
    
    m_isTrans = false;
}
//----- Printers -------------------------------------------------

void RedeBase::escreveX(const std::string &filePath, uint sampleSize, int *chosenNeurons, int header) const {
    uint i, n;
    bool allocedHere = false;

    if (xMat == NULL) {
    	Log::Warning("Não foi possível escrever x! Matriz não alocada!");
    	return;
    }

    std::ofstream os(filePath.c_str());

    if (os.fail()) {
        Log::Warning("Não foi possível abrir arquivo ", filePath.c_str());
        return;
    }

    if (header == 2)
    	os << stringDetalhes();

        //Se os neurons a serem escritos não foram especificados, escreve os primeiros
    if (chosenNeurons == NULL) {
        allocedHere = true;
        chosenNeurons = new int[sampleSize];

        for (i=0; i < sampleSize; i++)
            chosenNeurons[i] = i;
    }

        //Header
    if (header) {
        os << "step";
        for (i=0; i < sampleSize; i++) {
            os << "\tN" << chosenNeurons[i];
        }
        os << std::endl;
    }

        //Escreve no arquivo
    for (n=0; n < numStep; n++) {
	    os << n;
	    for (i=0; i < sampleSize; i++) {
	        os << '\t' << xMat[ chosenNeurons[i] ][n];
	    }
	    os << std::endl;
    }
    

    if (allocedHere)
        delete[] chosenNeurons;

    os.close();
}

void RedeBase::escreveY(const std::string &filePath, uint sampleSize, int *chosenNeurons, int header) const {
    uint i, n;
    bool allocedHere = false;

    if (yMat == NULL) {
    	Log::Warning("Não foi possível escrever y! Matriz não alocada!");
    	return;
    }

    std::ofstream os(filePath.c_str());

    if (os.fail()) {
        Log::Warning("Não foi possível abrir arquivo ", filePath.c_str());
        return;
    }

    if (header == 2)
    	os << stringDetalhes();

        //Se os neurons a serem escritos não foram especificados, escreve os primeiros
    if (chosenNeurons == NULL) {
        allocedHere = true;
        chosenNeurons = new int[sampleSize];

        for (i=0; i < sampleSize; i++)
            chosenNeurons[i] = i;
    }

        //Header
    if (header) {
        os << "step";
        for (i=0; i < sampleSize; i++) {
            os << "\tN" << chosenNeurons[i];
        }
        os << std::endl;
    }

        //Escreve no arquivo
    for (n=0; n < numStep; n++) {
	    os << n;
	    for (i=0; i < sampleSize; i++) {
	        os << '\t' << yMat[ chosenNeurons[i] ][n];
	    }
	    os << std::endl;
    }


    if (allocedHere)
        delete[] chosenNeurons;

    os.close();
}

void RedeBase::escreveXY(const std::string &filePath, uint sampleSize, int *chosenNeurons, int header) const {
    uint i, n;
    bool allocedHere = false;

    //----- Checagem de Erros -----------------------------------
    if (xMat == NULL || yMat == NULL) {
    	Log::Warning("Não foi possível escrever xy! xMat ou yMat não alocados!");
    	return;
    }

    std::ofstream os(filePath.c_str());

    if (os.fail()) {
        Log::Warning("Não foi possível abrir arquivo ", filePath.c_str());
        return;
    }
    //-----------------------------------------------------------


    if (header == 2)
    	os << stringDetalhes();


        //Se os neurons a serem escritos não foram especificados, escreve os primeiros
    if (chosenNeurons == NULL) {
        allocedHere = true;
        chosenNeurons = new int[sampleSize];

        for (i=0; i < sampleSize; i++)
            chosenNeurons[i] = i;
    }


        //Header
    if (header) {
        os << "step";
        for (i=0; i < sampleSize; i++) {
            os << "\tx_N" << chosenNeurons[i] << "\ty_N" << chosenNeurons[i];
        }
        os << std::endl;
    }

        //Escreve no arquivo
    for (n=0; n < numStep; n++) {
	    os << n;
	    for (i=0; i < sampleSize; i++) {
	        os << '\t' << xMat[ chosenNeurons[i] ][n] << '\t' << yMat[ chosenNeurons[i] ][n];
	    }
	    os << std::endl;
    }


    if (allocedHere)
        delete[] chosenNeurons;

    os.close();
}

void RedeBase::escreveFase(const std::string &filePath, uint sampleSize, int *chosenNeurons, int header) const {
    uint i, n;
    bool allocedHere = false;

    if (fase == NULL) {
        Log::Warning("Não foi possível escrever fase. Vetor não alocado!");
        return;
    }

    std::ofstream outFile(filePath.c_str());

    if (outFile.fail()) {
        Log::Warning("Não foi possível abrir arquivo ", filePath.c_str());
        return;
    }

    if (header == 2)
    	outFile << stringDetalhes();

        //Se os neurons a serem escritos não foram especificados, escreve os primeiros
    if (chosenNeurons == NULL) {
        allocedHere = true;
        chosenNeurons = new int[sampleSize];

        for (i=0; i < sampleSize; i++)
            chosenNeurons[i] = i;
    }

        //Header
    if (header) {
        outFile << "step";
        for (i=0; i < sampleSize; i++) {
            outFile << '\t' << "N" << chosenNeurons[i];
        }
        outFile << std::endl;
    }

        //Para cada step dentro do intervalo válido
    for (n=faseMin; n < faseMax; n++) {
        outFile << n - faseMin;

            //Para cada neuronio, escreve fase do step
        for (i=0; i < sampleSize; i++) {
            outFile << '\t' << fase[ chosenNeurons[i] ][n];
        }

        outFile << std::endl;
    }

    outFile.close();   

    if (allocedHere)
        delete[] chosenNeurons;
}

void RedeBase::escreveBurstStart(const std::string &filePath, uint sampleSize, int *chosenNeurons, int header) const {
    uint i, k;
    bool allocedHere = false;

    if (burstStart == NULL) {
        Log::Warning("Não é possível escrever burstStart. Vetor não alocado!");
        return;
    }

    std::ofstream outFile(filePath.c_str());

    if (outFile.fail()) {
        Log::Warning("Não foi possível abrir arquivo ", filePath.c_str());
        return;
    }

    if (header == 2) {
        outFile << stringDetalhes();
        outFile << "#Tempos nos quais ocorreram inícios de burst.\n";
    }

        //Se os neurons a serem escritos não foram especificados, escreve os primeiros
    if (chosenNeurons == NULL) {
        allocedHere = true;
        chosenNeurons = new int[sampleSize];

        for (i=0; i < sampleSize; i++)
            chosenNeurons[i] = i;
    }

        //Header
    if (header) {
        for (i=0; i < sampleSize; i++)
            outFile << "N" << chosenNeurons[i] << '\t';
        
        outFile << std::endl;
    }


    bool noneLeft = false;
    k=0;

    while (!noneLeft) {
        noneLeft = true;

        for (i=0; i < sampleSize; i++) {

            if (k < burstStart[i].size()) {
                outFile << burstStart[ chosenNeurons[i] ].at(k) << '\t';
                noneLeft = false;
            }
            else {
                outFile << 0 << '\t';
            }

        }

        outFile << std::endl;
        k++;
    }

    outFile.close();   

    if (allocedHere)
        delete[] chosenNeurons;
}

void RedeBase::escreveOrdem(const std::string &filePath, bool writeStep, int header) const {
    //----- Checagem de Erro -----
    if (ordem == NULL) {
        Log::Warning("Não foi possível escrever ordem. Vetor não alocado!");
        return;
    }
    
    if (faseMax <= faseMin) {
        Log::Warning("Não foi possível escrever ordem! O valor máximo da fase não é maior que o mínimo!");
        return;
    }
    //----------------------------
    

    std::ofstream outFile(filePath.c_str());

    if (outFile.fail()) {
        Log::Warning("Não foi possível abrir arquivo ", filePath.c_str());
        return;
    }

    if (header == 2)
    	outFile << stringDetalhes();

    //Imprime header
    if (header) {
        if (writeStep)
            outFile << "step\tordem\n";
        else
            outFile << "ordem\n";
    }
    
    if (writeStep) {
        for (uint n=faseMin; n < faseMax; n++)
            outFile << n << '\t' << ordem[n] << std::endl;;
    }
    else {
        for (uint n=faseMin; n < faseMax; n++)
            outFile << ordem[n] << std::endl;;   
    }
    
    outFile.close();
}

void RedeBase::escreveCampoMed(const std::string &filePath, bool writeStep, int header) const {
	if (campoMed == NULL) {
        Log::Warning("Não é possível escrever campoMed. Vetor não alocado!");
        return;
    }

    std::ofstream outFile(filePath.c_str());

    if (outFile.fail()) {
        Log::Warning("Não foi possível abrir arquivo ", filePath.c_str());
        return;
    }

    if (header == 2)
    	outFile << stringDetalhes();

    //Imprime header
    if (header) {
        if (writeStep)
            outFile << "step\tcampoMed\n";
        else
            outFile << "campoMed\n";
    }
    
    if (writeStep) {
        for (uint n=0; n < numStep; n++)
            outFile << n << '\t' << campoMed[n] << std::endl;;
    }
    else  {
        for (uint n=0; n < numStep; n++)
            outFile << campoMed[n] << std::endl;;   
    }
    
    outFile.close();
}

//----- Exporter and Imporer ----------------------------------------
void RedeBase::exportRede(const std::string &filePath) const{
    uint i, n;
    
    std::ofstream os(filePath.c_str());

    //----- Dados Estáticos ------------------------------------------------------
    os << "#Tipo da Rede\n";
    os << "Global\n";

    os << "#Numero de Neuronios\n";
    os << numNeurons << std::endl;

    os << "#Numero de Steps\n";
    os << numStep << std::endl;

    os << "#Transiente\n";
    os << transiente << std::endl;

    os << "#Eps\n";
    os << eps << std::endl;

    os << "#Beta\n";
    os << beta << std::endl;

    os << "#Sigma\n";
    os << sig << std::endl;

    os << "#faseMin\n";
    os << faseMin << std::endl;

    os << "#faseMax\n";
    os << faseMax << std::endl;

    os << "#Parâmetro de Ordem Médio\n";
    os << ordemMed << std::endl;

    os << "#Média do Campo Médio\n";
    os << avgCampoMed << std::endl;

    os << "#Variância do Campo Médio\n";
    os << varCampoMed << std::endl;

    os << "#m_shouldCalc\n";
    for (i=0; i < CALC_LIST_SIZE; i++)
         os << m_shouldCalc[i] << std::endl;

    if (m_alpha == NULL) {
        std::cout << "Can't save alpha. Vector not allocated!\n";
        exit(1);
    }
    os << "#Alpha de cada Neurônio\n";
    for (i=0; i < numNeurons; i++)
        os << m_alpha[i] << ' ';
    os << std::endl;
    //----------------------------------------------------------------------------


        //xMat
    if (m_shouldCalc[Grandeza::x]) {
        if (xMat == NULL) {
            std::cout << "Can't save xMat. Matrix not allocated!\n";
            exit(1);
        }

        os << "#xMat  -  Matriz com x para cada Neuronio para cada tempo.\n";

        for (i=0; i < numNeurons; i++) {
            for (n=0; n < numStep; n++) {
                os << xMat[i][n] << ' ';
            }
            os << std::endl;
        }
    }

        //yMat
    if (m_shouldCalc[Grandeza::y]) {
        if (yMat == NULL) {
            std::cout << "Can't save yMat. Matrix not allocated!\n";
            exit(1);
        }

        os << "#yMat  -  Matriz com y para cada Neuronio para cada tempo.\n";

        for (i=0; i < numNeurons; i++) {
            for (n=0; n < numStep; n++) {
                os << yMat[i][n] << ' ';
            }
            os << std::endl;
        }
    }

        //Fase
    if (m_shouldCalc[Grandeza::fase]) {
        if (fase == NULL) {
            std::cout << "Can't save fase. Matrix not allocated!\n";
            exit(1);
        }

        os << "#fase  -  Matrix com fase para cada Neuronio para cada tempo.\n";

        for (i=0; i < numNeurons; i++) {
            for (n=0; n < numStep; n++) {
                os << fase[i][n] << ' ';
            }
            os << std::endl;
        }
    }

        //Ordem
    if (m_shouldCalc[Grandeza::ordem]) {
        if (ordem == NULL) {
            std::cout << "Can't save ordem. Vector not allocated!\n";
            exit(1);
        }

        os << "#ordem  -  Vetor com Parâmetro de Ordem para cada tempo.\n";

        for (n=0; n < numStep; n++)
            os << ordem[n] << ' ';
        os << std::endl;
    }

        //Campo Médio
    if (m_shouldCalc[Grandeza::campoMed]) {
        if (campoMed == NULL) {
            std::cout << "Can't save campoMed. Vector not allocated!\n";
            exit(1);
        }

        os << "#campoMed  -  Vetor com o Campo Médio para cada tempo.\n";

        for (n=0; n < numStep; n++)
            os << campoMed[n] << ' ';
        os << std::endl;
    }


    os.close();
}

void RedeBase::importRede(const std::string &filePath) {
    uint i, n;
    std::string line;
    
    std::ifstream is(filePath.c_str());

    //----- Checagem de Erro -----
    if (is.fail()) {
        Log::Error("Não foi possível importar rede do arquivo ", filePath.c_str());
        exit(1);
    }

    std::getline(is, line); std::getline(is, line);
        //Checa tipo de rede
    if (line != "Global") {
        Log::Error("Failed to load Rede. Expected tipo Global, but saved Rede is ", line.c_str());
        exit(1);
    }
    //----------------------------

        //Limpa rede atual
    destroy();

    //----- Dados Estáticos ------------------------------------------------------
        //numNeurons
    std::getline(is, line);
    is >> numNeurons;
    is.ignore(1);

        //numStep
    std::getline(is, line);
    is >> numStep;
    is.ignore(1);

        //Transiente
    std::getline(is, line);
    is >> transiente;
    is.ignore(1);

        //Eps
    std::getline(is, line);
    is >> eps;
    is.ignore(1);

        //Beta
    std::getline(is, line);
    is >> beta;
    is.ignore(1);

        //Sigma
    std::getline(is, line);
    is >> sig;
    is.ignore(1);

        //faseMin
    std::getline(is, line);;
    is >> faseMin;
    is.ignore(1);

        //faseMax
    std::getline(is, line);
    is >> faseMax;
    is.ignore(1);

        //ordemMed
    std::getline(is, line);
    is >> ordemMed;
    is.ignore(1);

        //avgCampoMed
    std::getline(is, line);
    is >> avgCampoMed;
    is.ignore(1);

        //varCampoMed
    std::getline(is, line);
    is >> varCampoMed;
    is.ignore(1);


        //m_shouldCalc
    std::getline(is, line);
    for (i=0; i < CALC_LIST_SIZE; i++) {
         is >> m_shouldCalc[i];
         is.ignore(1);
    }


        //Alpha
    m_alpha = new double[numNeurons];
    std::getline(is, line);
    for (i=0; i < numNeurons; i++)
        is >> m_alpha[i];
    is.ignore(2);
    //----------------------------------------------------------------------------


        //xMat
    if (m_shouldCalc[Grandeza::x]) {
        xMat = allocMat<float>(numNeurons, numStep);

        std::getline(is, line);

        for (i=0; i < numNeurons; i++) {
            for (n=0; n < numStep; n++) {
                is >> xMat[i][n];
            }       
        }
        is.ignore(2);
    }

        //yMat
    if (m_shouldCalc[Grandeza::y]) {
        yMat = allocMat<float>(numNeurons, numStep);

        std::getline(is, line);

        for (i=0; i < numNeurons; i++) {
            for (n=0; n < numStep; n++) {
                is >> xMat[i][n];
            }       
        }
        
        is.ignore(2);
    }

        //Fase
    if (m_shouldCalc[Grandeza::fase]) {
        fase = allocMat<float>(numNeurons, numStep);

        std::getline(is, line);

        for (i=0; i < numNeurons; i++) {
            for (n=0; n < numStep; n++) {
                is >> fase[i][n];
            }       
        }
        is.ignore(2);
    }

        //Ordem
    if (m_shouldCalc[Grandeza::ordem]) {
        ordem = new double[numStep];

        std::getline(is, line);

        for (n=0; n < numStep; n++)
            is >> ordem[n];
        
        is.ignore(2);
    }

        //Campo Médio
    if (m_shouldCalc[Grandeza::campoMed]) {
        campoMed = new double[numStep];

        std::getline(is, line);

        for (n=0; n < numStep; n++)
            is >> campoMed[n];
    }


    is.close();
}


//----- Setters --------------------------------------------------------


void RedeBase::setSynParams(double limiar, double potencial) {
    limiarSinapse = limiar;
    potencialSinapse = potencial;
}

void RedeBase::setCalc(enum Grandeza::Grandeza grandEnum, bool shouldCalc) {
    //Setta grandeza pedida
    m_shouldCalc[grandEnum] = shouldCalc;

        //Se o usuario mandou Grandeza::all, setta todos
    if (grandEnum == Grandeza::all) {
        for (uint i=0; i < CALC_LIST_SIZE; i++)
            m_shouldCalc[i] = shouldCalc;

        m_shouldCalc[Grandeza::onlyOrdemMed]    = false;
        m_shouldCalc[Grandeza::onlyAvgCampoMed] = false;
    }

        //Se o pedido é acionar onlyOrdemMed, aciona também ordem
    if (grandEnum == Grandeza::onlyOrdemMed && shouldCalc == true)
        m_shouldCalc[Grandeza::ordem] = true;

        //Se o pedido é acionar onlyAvgCampoMed, aciona também campoMed
    if (grandEnum == Grandeza::onlyAvgCampoMed && shouldCalc == true)
        m_shouldCalc[Grandeza::campoMed] = true;

    	//Se o pedido é a variância do campoMed, precisa calcular o campoMed
    if (grandEnum == Grandeza::varCampoMed && shouldCalc == true)
    	m_shouldCalc[Grandeza::campoMed] = true;


        //calcula() não consegue mandar calcular só o y
    if (grandEnum == Grandeza::y)
        m_shouldCalc[Grandeza::x] = true;
}

void RedeBase::ligarSinal() { 
    if (!m_sinal->isSet()) {
        Log::Error("Sinal externo ligado sem ser settado!");
        exit(1);
    }

    m_hasExternalSignal = true;
}


//----- Allocators -----------------------------------------------------

void RedeBase::allocForFase() {
    if (m_x == NULL)
        m_x = new double[numNeurons];

    if (m_y == NULL)
        m_y = new double[numNeurons];

    if (m_burstId == NULL)
        m_burstId = new uint[numNeurons];
        
    if (m_histX == NULL)
    	m_histX = allocMat<float>(numNeurons, HIST_X_SIZE);
}

void RedeBase::allocForOrdem() {
    if (m_real == NULL)
        m_real = new double[numStep];
        
    if (m_imag == NULL)
        m_imag = new double[numStep];
}

//----- Libera memória -------------------------------------------------

void RedeBase::freeFase() {
    if (fase == NULL)  return;

        //Para cada neurônio, libera a memória de fase
    for (uint i=0; i < numNeurons; i++) {
        delete[] fase[i];
    }

    delete[] fase;

    fase = NULL;
}

void RedeBase::freeOrdem() {
    if (ordem != NULL) {
        delete[] ordem;
        ordem = NULL;
    }
}

void RedeBase::freeCampoMed() {
    if (campoMed != NULL) {
        delete[] campoMed;
        campoMed = NULL;
    }
}

void RedeBase::freeXMat() {
	if (xMat != NULL) {
		for (uint i=0; i < numNeurons; i++)
			delete[] xMat[i];
		delete[] xMat;

		xMat = NULL;
	}
}

void RedeBase::freeYMat() {
	if (yMat != NULL) {
		for (uint i=0; i < numNeurons; i++)
			delete[] yMat[i];
		delete[] yMat;

		yMat = NULL;
	}
}

void RedeBase::freeXYMat() {
    freeXMat();
    freeYMat();
}

void RedeBase::freeBurstStart() {
    if (burstStart != NULL)
        delete[] burstStart;

    burstStart = NULL;
}

void RedeBase::freeAuxiliares() {
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
}

void RedeBase::destroy() {
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

    clearCalcList();
}


//----- Auxiliares ------------------------------------------------------


void RedeBase::setInitialCondition() {
    uint i;
    srand(ciSeed);
    
    size_t x_size = m_initCond[0].size();
    size_t y_size = m_initCond[1].size();
    
    
    //----- x ----------
    //Random
    if (x_size == 0) {
        for (i=0; i < numNeurons; i++)
            m_x[i] = -2.0  +  4.0 * (rand() % 10000000)/10000000.0;
    }
    
    //Caso em que há um x0 estabelecido para cada neurônio
    else if (x_size == numNeurons) {
        for (i=0; i < numNeurons; i++)
            m_x[i] = m_initCond[0].at(i);
    }
    
    //Caso em que apenas alguns valores de x0 foram estabelecidos
    else {
        for (i=0; i < numNeurons; i++)
            m_x[i] = m_initCond[0].at( i % x_size );
    }
    //------------------
    
    
    //----- y ----------
    //Caso Random
    if (y_size == 0) {
        for (i=0; i < numNeurons; i++)
            m_y[i] = -4.0  +  4.0 * (rand() % 10000000)/10000000.0;
    }

    //Caso em que há um y0 estabelecido para cada neurônio
    else if (y_size == numNeurons) {
        for (i=0; i < numNeurons; i++)
            m_y[i] = m_initCond[1].at(i);
    }

    //Caso em que apenas alguns valores de y0 foram estabelecidos
    else {
        for (i=0; i < numNeurons; i++)
            m_y[i] = m_initCond[1].at( i % x_size );
    }
    //------------------
}

bool RedeBase::checaBurst(uint i) {
    //Checa se o tempo atual é o inicio de um burst

        //Se x for menor que o limiar, não é inicio de burst    
    if (m_x[i] < LIMIAR_SPIKE)
        return false;
    
        //Se algum tempo anterior próximo foi burst, não é inicio de burst    
    for (int k=0; k < HIST_X_SIZE; k++) {
        if (m_histX[i][k] > LIMIAR_SPIKE)
            return false;
    }

    return true;
}

const std::string RedeBase::stringDetalhes() const {
	std::string str = "#NumNeurons = " + toString(numNeurons) + "  --  NumStep = " + toString(numStep) + "  --  Transiente = " + toString(transiente) + "\n#Alpha = ";
	

	//Alpha
	if (m_alpha[1] != m_alpha[2] && m_alpha[1] != m_alpha[3]) {
		str = str + "Rand";
	}
	else {
		str = str + toString(m_alpha[1]);
	}

	//----- Condições Iniciais -----
	str = str + "  --  (x0,y0) = (";
	if (m_initCond[0].size() == 0) {
		str = str + "Rand, ";
	}
    else if (m_initCond[0].size() == 1) {
        str = str + toString(m_initCond[0][0]) + ", ";
    }
	else {
		str = str + "set-by-hand, ";
	}

    if (m_initCond[1].size() == 0) {
        str = str + "Rand)\n";
    }
    else if (m_initCond[1].size() == 1) {
        str = str + toString(m_initCond[1][0]) + ")\n";
    }
    else {
        str = str + "set-by-hand)\n";
    }
    //------------------------------

    return str;
}
