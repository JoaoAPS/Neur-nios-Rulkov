#ifndef REDE_H
#define REDE_H

#include "Base.h"


class Rede : public RedeBase {
public:
	uint adjVetSize;
    double conectividadeMedia;
    uint *adjVet;

protected:
		//Guarda a contribuição dos vizinhos à dinâmica de cada neurônio, para um determinado tempo
	double *m_contribVizinhos;

public:
    Rede(uint numNeurons, uint numStep, uint transiente, double eps=0.0, double beta=0.001, double sig=0.001);
    Rede(uint numNeurons, uint numStep, uint transiente, double eps, const std::string &adjVetPath);
    virtual ~Rede();

    //Calculators
    virtual void calcula();

	//Setters
	void readAdjVet(const std::string &filePath);
    void setAdjVet(uint *_adjVet, int _adjVetSize, double _conectividadeMedia);

    //Printers
    void escreveAdjMat(const std::string &filePath, bool asPGM=false);

	//Memória
    void freeAuxiliares();
    virtual void destroy();

private:
    virtual void evoluiStep(uint n);
    virtual void calcContribVizinhos();
    virtual void calcContribVizinhos_eletrico();
};


#endif
