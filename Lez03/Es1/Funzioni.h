#ifndef __funzioni_h__
#define __funzioni_h__

#include <iostream>
#include <cmath>
#include <vector>

//SE da problemi in fase di copmilazione commenta le funzioni print

/*
	DICHIARAZIONE FUNZIONE BASE - FUNZIONE MADRE ASTRATTA
*/

using namespace std;

class FunzioneBase{

public:

  virtual double Eval(double x) const = 0;
};


/*=======================================================================================
	DICHIARAZIONE FUNZIONE CALL-OPTION
========================================================================================*/

class CallOption : public FunzioneBase {

public:

	CallOption(double T, double sigma, double K, double St, double r);

	~CallOption() {;};

  double GetT() const {return m_T;};
  double GetSigma() const {return m_sigma;};
  double GetK() const {return m_K;};
  double GetSt() const {return m_St;};
  double Getr() const {return m_r;};

	virtual double Eval(double t) const;
	
private:
	double m_T, m_sigma, m_K, m_St, m_r;

};

/*=======================================================================================
	DICHIARAZIONE FUNZIONE PUT-OPTION
========================================================================================*/

class PutOption : public FunzioneBase {

public:

	PutOption(double T, double sigma, double K, double St, double r);

	~PutOption() {;};

  double GetT() const {return m_T;};
  double GetSigma() const {return m_sigma;};
  double GetK() const {return m_K;};
  double GetSt() const {return m_St;};
  double Getr() const {return m_r;};

	virtual double Eval(double t) const;
	
private:
	double m_T, m_sigma, m_K, m_St, m_r;

};

#endif
