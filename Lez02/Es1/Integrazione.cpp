#include "Integrazione.h"
#include "Funzioni.h"
#include "random.h"
//#include "statistica.h"
#include <cmath>
#include <vector>

using namespace std;

Integral::~Integral() {};

//=============================================================================

//Metodi MEDIA
double Integral::Media(const FunzioneBase* f, double xmin, double xmax, int punti, Random& m_gen) {

	double sum = 0;
	double x = 0;
	
	for (int i=0; i<punti; i++) {
		x = m_gen.Rannyu(xmin, xmax);
		double fi = f->Eval(x);
		sum = sum + fi;
	}

	return (sum/punti)*(xmax-xmin);
};


/*========================================================================================
IMPORTANCE SAMPLING
=========================================================================================*/
double Integral::ImpSampling(FunzioneBase* f, FunzioneBase* pdf, double xmin, double xmax, int punti, Random& m_gen) {
	
	double sum=0;
	double x = 0;
	
	for (int i=0; i<punti; i++){
		x = m_gen.Campionamento(xmin, xmax);
		double fi = (f->Eval(x))/(pdf->Eval(x));
		sum = sum + fi;
	}
	
	return (sum/punti)*(xmax-xmin);
};

