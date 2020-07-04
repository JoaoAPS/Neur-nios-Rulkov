#ifndef GLOBAL_H
#define GLOBAL_H

#include "Base.h"


class RedeGlobal : public RedeBase {
public:     //Métodos
    RedeGlobal(uint numNeurons, uint numStep, uint transiente, float eps=0.0, float beta=0.001, float sig=0.001);
    ~RedeGlobal();

protected:
    void evoluiStep(uint n);
    void evoluiStep_eletrico(uint n);
};

#endif
