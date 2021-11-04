/// closest octave, only kr

#include "SC_PlugIn.h"

static InterfaceTable *ft;

struct Nextoct : public Unit {
  
};


extern "C"
{
void load(InterfaceTable *inTable);
void Nextoct_next_k(Nextoct *unit, int inNumSamples);
void Nextoct_Ctor(Nextoct* unit);
};


void Nextoct_Ctor(Nextoct* unit)
{
  SETCALC(Nextoct_next_k);
  Nextoct_next_k(unit, 1);
}

// inCandidate has to smaller than inFreq !!

void Nextoct_next_k(Nextoct *unit, int inNumSamples) {
  float *out = OUT(0);
  float inFreq = IN0(0);
  float inCandidate = IN0(1);
  
  float candidate = inFreq;
  float lastCandidate = inFreq;

  int harmCounter = 1;
  
  for (int i=0; i < inNumSamples; ++i) {
    while (sc_abs(candidate - inCandidate) <= sc_abs(lastCandidate - inCandidate)) {
      lastCandidate = candidate;
      if (inCandidate > inFreq) {
	candidate = inFreq * harmCounter;
      } else {
	candidate = inFreq / harmCounter;
      }
      harmCounter = harmCounter * 2;
    }
    out[i] = lastCandidate;
  }
}


PluginLoad(Nextoct)
{
    // InterfaceTable *inTable implicitly given as argument to the load function
    ft = inTable; // store pointer to InterfaceTable

    DefineSimpleUnit(Nextoct);
}
