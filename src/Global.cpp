#include <include/Global.h>


//----- Constructors ------------------------------------------

RedeGlobal::RedeGlobal(uint numNeurons, uint numStep, uint transiente, float eps, float beta, float sig) :
RedeBase(numNeurons, numStep, transiente, eps, beta, sig) { }

RedeGlobal::~RedeGlobal() {
    destroy();
}


//----- Calculators -------------------------------------------

void RedeGlobal::evoluiStep_eletrico(uint n) {
    uint i;
    double media, tmp;

    //Calcula a média de x
    media = 0;
    for (i=0; i < numNeurons; i++) {
        media += m_x[i];
    }
    media /= numNeurons;

    //Calcula x e y
    for (i=0; i < numNeurons; i++) {
        tmp = m_x[i];

        m_x[i] = m_alpha[i] / (1.0 + m_x[i] * m_x[i]) + m_y[i] + eps * media;
        m_y[i] = m_y[i] - sig * tmp - beta;
    }

    //Aplica o sinal externo se houver    
    if (m_hasExternalSignal)
        if ((int)n >= m_inicioSinal)
            m_sinal->aplica(m_x, m_y, n);


    //Chama função externa do usuário
    externalFunction(n);
}


void RedeGlobal::evoluiStep(uint n) {
    uint i, j;
    double tmp;
    
    //Calcula efeito das sinapses
    for (i=0; i < numNeurons; ++i) {
        //Se o neurônio atual não está disparando, não influencia os outros
        if (m_x[i] < limiarSinapse)
            continue;
        
        //Afeta cada outro neurônio
        for (j=0; j < numNeurons; j++) {
            if (i == j)
                continue;
            
            m_x[j] += (eps / (numNeurons-1)) * (potencialSinapse - m_x[j]);
        }
    }

    //Calcula x e y
    for (i=0; i < numNeurons; ++i) {
        tmp = m_x[i];

        m_x[i] = m_alpha[i] / (1.0 + m_x[i] * m_x[i]) + m_y[i];
        m_y[i] = m_y[i] - sig * tmp - beta;
    }

    //Aplica o sinal externo se houver    
    if (m_hasExternalSignal)
        if ((int)n >= m_inicioSinal)
            m_sinal->aplica(m_x, m_y, n);


    //Chama função externa do usuário
    externalFunction(n);
}
