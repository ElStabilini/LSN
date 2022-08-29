#ifndef __Integrazione_h__
#define __Integrazione_h__

#include "random.h"
#include "Funzioni.h"
#include <cmath>
#include <vector>

using namespace std;

class Integral {
	
public :

//	Integral(unsigned int seed);
	~Integral();

	double Media(const FunzioneBase* f, double xmin, double xmax, int punti, Random& m_gen);
	double ImpSampling(FunzioneBase* f, FunzioneBase* pdf, double xmin, double xmax, int punti, Random& m_gen);
	
//IMPORTANCE SAMPLING

	
	double GetErrore() const {return m_errore;}
	unsigned int GetN() const {return m_punti;}

private :

	double m_errore;
	unsigned int m_punti;
};

#endif //__INTEGRALEMC_H__	
