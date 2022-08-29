#include "Funzioni.h"
#include <iostream>
#include <cmath>

/*=============================================================
IMPLEMENTAZIONE METODI FUNZIONE CALL-OPTION
=============================================================*/
CallOption::CallOption(double T, double sigma, double K, double St, double r) {
	m_T = T;
	m_sigma = sigma;
	m_K = K;
	m_St = St;
	m_r = r;
};


double CallOption::Eval(double t) const {

	double Ct = 0;
	if (t>=0 && t <=1) {
		double d1 = (log(m_St/m_K) + (m_r + 0.5*(m_T-t)*pow(m_sigma,2)) )/(m_sigma*sqrt(m_T-t));
		double d2 = d1 - m_sigma*sqrt(m_T - t);
		double N1 = 0.5*(1 + erf(d1/sqrt(2)) );
		double N2 = 0.5*(1 + erf(d2/sqrt(2)) );
		Ct = m_St*N1 - m_K*exp(-(m_T-t)*m_r)*N2;
	} else 
		cout << "modello non adatto" << endl;
		
	return Ct;
};

/*=============================================================
IMPLEMENTAZIONE METODI FUNZIONE PUT-OPTION
=============================================================*/

PutOption::PutOption(double T, double sigma, double K, double St, double r) {
	m_T = T;
	m_sigma = sigma;
	m_K = K;
	m_St = St;
	m_r = r;
};

double PutOption::Eval(double t) const {
	double Pt = 0;
	if (t>=0 && t <=1) {
		double d1 = (log(m_St/m_K) + (m_r + 0.5*(m_T-t)*pow(m_sigma,2)) )/(m_sigma*sqrt(m_T-t));
		double d2 = d1 - m_sigma*sqrt(m_T - t);
		double N1 = 0.5*(1 + erf(d1/sqrt(2)) );
		double N2 = 0.5*(1 + erf(d2/sqrt(2)) );
		Pt = m_St*(N1-1) - m_K*exp(-(m_T-t)*m_r)*(N2-1);
	} else 
		cout << "modello non adatto" << endl;
		
	return Pt;
}
