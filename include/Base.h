#ifndef BASE_H
#define BASE_H

#include "../../pch/preCompiledHeader.h"
#include "SinalExterno.h"


	//Parâmetros para encontrar fase
#define LIMIAR_SPIKE 0.0    //Quão grande deve ser a variável x para que um neurônio esteja disparando
#define HIST_X_SIZE 50      //Quanto tempo o neurônio deve ficar em silêncio antes de começar um novo burst
	//Ajustar conforme mais grandezas podem ser calculadas
#define CALC_LIST_SIZE 12


namespace Grandeza { enum Grandeza { none, all ,x, y, fase, ordem, campoMed, varCampoMed, burstStart, onlyOrdemMed, onlyAvgCampoMed }; }


class RedeBase {

//Variáveis
public:
    //Parâmetros Base
    uint numNeurons, numStep, transiente;
    double eps, beta, sig;
    int alphaSeed, ciSeed;
    //Características da Rede calculadas
    uint faseMin, faseMax;
    double *ordem, *campoMed;
    double ordemMed, avgCampoMed, varCampoMed;
    //Características do Neurônio
    float **xMat, **yMat, **fase;
    std::vector<uint> *burstStart;
    //Parâmetros da sinapse química
    double limiarSinapse, potencialSinapse;



protected:
    std::vector< std::vector<double> > m_initCond;
    double *m_alpha;
        //Guardam quais grandezas devem ser salvas em memória
    bool m_shouldCalc[CALC_LIST_SIZE];
    bool m_cleanUp;
        //Variáveis auxiliares para calculo de fase e ordem
    double *m_x, *m_y;
    uint  *m_burstId;
    float **m_histX;
    double *m_real, *m_imag;
    	//Variáveis do Sinal Externo
    SinalExterno* m_sinal;
    bool m_hasExternalSignal;
    int m_inicioSinal;
        //Auxiliar
    bool m_isTrans;

//Métodos
public: 
    RedeBase(uint numNeurons, uint numStep, uint transiente, double eps=0.0, double beta=0.001, double sig=0.001);
    virtual ~RedeBase();

        //Inits
    void init(std::vector<double> alpha=std::vector<double>(0), std::vector<double> x0=std::vector<double>(0), std::vector<double> y0=std::vector<double>(0));
    void init(float alpha, std::vector<double> x0=std::vector<double>(0), std::vector<double> y0=std::vector<double>(0));
    void init(std::vector<double> alpha, float x0, float y0);
    void init(float alpha, float x0, float y0); //init() chamado por float
    void clearCalcList();


	//Calculators
    virtual void calcula();
    virtual void externalFunction(uint n) { }


        //Printers
    void escreveX   (const std::string &filePath, uint sampleSize, int *chosenNeurons=NULL, int header=1)       const;
    void escreveY   (const std::string &filePath, uint sampleSize, int *chosenNeurons=NULL, int header=1)       const;
    void escreveXY  (const std::string &filePath, uint sampleSize, int *chosenNeurons=NULL, int header=1)       const;
    void escreveFase(const std::string &filePath, uint sampleSize, int *chosenNeurons=NULL, int header=1)       const;
    void escreveBurstStart(const std::string &filePath, uint sampleSize, int *chosenNeurons=NULL, int header=1) const;
    void escreveOrdem   (const std::string &filePath, bool writeStep=true,  int header=0) const;
    void escreveCampoMed(const std::string &filePath, bool writeStep=false, int header=1) const;




        //Setters
    void setSynParams(double limiar, double potencial);
    void setEps(float eps) { this->eps  = eps;  }
    void setSeeds(int _alphaSeed, int _ciSeed) { alphaSeed = _alphaSeed; ciSeed = _ciSeed; }
    void setSeed(int seed) { setSeeds(seed, seed); }
    void setCalc(enum Grandeza::Grandeza grandEnum, bool shouldCalc=true);
    void setRulkovParams(double _beta, double _sigma) { beta = _beta; sig = _sigma; }
    void setCleanUp(bool shouldCleanUp=true) { m_cleanUp = shouldCleanUp; }
    //Sinal Externo
    void setSinal(SinalExterno* sinal) { m_sinal = sinal; ligarSinal(); }
    void ligarSinal();
    void desligarSinal() { m_hasExternalSignal = false; }


        //Getters
    bool getShouldCalc(enum Grandeza::Grandeza grandEnum) const { return m_shouldCalc[grandEnum]; }
    double getAlpha(uint index) const { return m_alpha[index]; }
    double* getAuxX() { return m_x; }
    double* getAuxY() { return m_y; }
    double  getAuxX(uint i) { return m_x[i]; }
    double  getAuxY(uint i) { return m_y[i]; }
    std::vector< std::vector<double> > getInitCond() { return m_initCond; }


        //Exporter and Importer
    void exportRede(const std::string &filePath) const;
    void importRede(const std::string &filePath);


        //Libera memória
    void freeFase();
    void freeOrdem();
    void freeCampoMed();
    void freeXMat();
    void freeYMat();
    void freeXYMat();
    void freeBurstStart();
    virtual void freeAuxiliares();
    virtual void destroy();


protected:
		//Calculators
    void calcNone();
    void calcXY();
    void calcBurstStart();
    void calcFase();
    void calcOrdem();
    void calcCampoMed();
    void calcFase_CampoMed();
    void calcOrdem_CampoMed();
    void calcBurstStart_CampoMed();
    void calcFaseDadoX();
    void calcOrdemDadoFase();
    void calcOrdemDadoX();
    void calcBurstStartDadoX();
    void calcCampoMedDadoX();
    void calcVarCampoMed();

    void calcTransiente();
    virtual void evoluiStep(uint n) = 0;


        //Allocators
    void allocForFase();
    void allocForOrdem();

        //Auxiliares
    void setInitialCondition();
    bool checaBurst(uint i);
    const std::string stringDetalhes() const;
};

#endif
