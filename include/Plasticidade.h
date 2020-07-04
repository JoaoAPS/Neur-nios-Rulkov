#ifndef PLASTICIDADE_H
#define PLASTICIDADE_H

#include "Rede.h"
#include <deque>

#define SYNWEIGHT_MIN 0.0   //Valor mínimo que a força de sinapse pode assumir (lembrando que a sinapse ainda será multiplicada por eps)
#define SYNWEIGHT_MAX 1.0   //Valor máximo que a força de sinapse pode assumir (lembrando que a sinapse ainda será multiplicada por eps)
#define NUM_SPIKES_PLAS 8   //Número de spikes mais recentes do neurônio i que influenciam um spike atual do neurônio j

namespace Modelo { enum Modelo {BTDP, Spike, Estatico, none}; }

class RedePlasticidade : public Rede {
public:
    //Variáveis da plasticidade
    std::vector<double> synWeight;
    double amplPotenc, amplDepres, tauPotenc, tauDepres;
    uint tempoSaturacao;
    

protected:
    enum Modelo::Modelo m_modelo;
    int *m_lastBurstStart;
    std::deque<uint> *m_lastSpikesTime;
    double m_amplPotInterna, m_amplDepInterna, m_angCoefBTDP;
    std::vector<uint> m_neuronsBurstingNow;
    bool attXHistOnPlasticity;


public:
    RedePlasticidade(uint numNeurons, uint numStep, uint transiente, double eps, enum Modelo::Modelo modelId=Modelo::none);
    ~RedePlasticidade();
    void destroy();

    //Calculators
    void calcula();
    void tiraCITrans(int duracao);
    
    //Setters
    void setPlasticityParams(double _amplPotenc, double _amplDepres, double _tauPotenc, double _tauDepres);
    void setPlasticityParams(double _amplPotenc, double _amplDepres, uint _tempoSaturacao);
    void setSynWeights(std::vector<double> _synWeight);
    void setSynWeights(double _synWeight) { setSynWeights(std::vector<double>(1, _synWeight)); }
    void readSynWeights(const std::string &weightFilePath);
    void setRulkovParamsBTDP() { beta  = 0.0011; sig = 0.0009; }

    //Writters
    void escreveSynWeights(const std::string &filePath, int header=1, bool writeAdjVet=true) const;
    

protected:
    void evoluiStep(uint n);
    void calcContribVizinhos();
    void atualizaSynWeight_BTDPModel(uint n);
    void atualizaSynWeight_spikeModel(uint n);
    double calcMudancaBTDP(int burstStartLatency);
};


#endif
