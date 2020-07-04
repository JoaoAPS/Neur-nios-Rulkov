#ifndef SINALHARMONICO_H
#define SINALHARMONICO_H

#include "../../pch/preCompiledHeader.h"


class SinalExterno {

typedef void(*SinalFuncPtr)(double&, double&, int, void*);

private:
	//Devem ser settados por setSinal()e setAfetados()
    uint m_numAfetados, *m_afetados;
    double m_fracAfetados;

    //Pointers to outside data
    void *m_params;			//Poiter to struct of parameters
    SinalFuncPtr m_func;	//Function pointer


public:
    uint startStep;


public:
    SinalExterno();
    ~SinalExterno();


    //Setters
    void setSinal(SinalFuncPtr func, void *params=NULL, int startStep=0);
    void setAfetados(uint numNeurons, double fracAfetados, int *chosenNeurons=NULL);


    //Calculators
    void aplica(double *x, double *y, int n);


	//Getters
    double getFracAfetados()   const { return m_fracAfetados; }
    uint getNumAfetados()      const { return m_numAfetados;  }
    uint getAfetado(int index) const { return m_afetados[index]; }
    bool isAfetadosSet()       const { if (m_afetados == NULL) return false; return true; }
    bool isOutsideDataSet()	   const { if (m_func != NULL) return true; return false; }
    bool isSet()			   const { if (isAfetadosSet() && isOutsideDataSet()) return true; return false; }
};


#endif