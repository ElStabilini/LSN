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
	virtual void Print() const = 0;
	//virtual double operator () (double x) const = 0;
};


/*============================================================================================
	DICHIARAZIONE FUNZIONE PARABOLA
============================================================================================*/
class Parabola: public FunzioneBase{

public:
	
	//costruttori
  Parabola();
  Parabola(double a, double b, double c) {m_a = a; m_b = b; m_c = c;};

	//distruttore
  ~Parabola() {;};

  void SetA(double a) {m_a = a;};
  void SetB(double b) {m_b = b;};
  void SetC(double c) {m_c = c;};
  double GetA() const {return m_a;};
  double GetB() const {return m_b;};
  double GetC() const {return m_c;};
  double Inversa(double x) const {return (1-sqrt(1-x));};

	virtual double Eval(double x) const {return m_a*x*x + m_b*x + m_c;};
	virtual void Print () const;


private:

  double m_a, m_b, m_c;

};



/* =========================================================================
	DICHIARAZIONE FUNZIONE COSENO
==========================================================================*/

class Coseno : public FunzioneBase {

public:
	Coseno() {m_a = 0; m_b = 1; m_c = 0; m_d = 0;};
	Coseno(double a, double b, double c, double d) { m_a = a;	m_b = b; m_c = c; m_d = d;};

  ~Coseno() {;};

	void SetA(double a) {m_a = a;};
  void SetB(double b) {m_b = b;};
  void SetC(double c) {m_c = c;};
	void SetD(double d) {m_d = d;};
  double GetA() const {return m_a;};
  double GetB() const {return m_b;};
 	double GetC() const {return m_c;};
	double GetD() const {return m_d;};

	virtual double Eval(double x) const {return m_a*cos(m_b*x + m_c) + m_d;};

	virtual void Print() const;

private:
	double m_a, m_b, m_c, m_d;

};


/*=======================================================================================
	DICHIARAZIONE FUNZIONE RETTA
========================================================================================*/

class Retta : public FunzioneBase {

public:

	Retta() {m_a = 0; m_b = 0;};
	Retta(double a, double b) {m_a = a; m_b = b;};

	~Retta() {;};

  double GetA() const {return m_a;};
  double GetB() const {return m_b;};
	void SetA(double a) {m_a = a;};
  void SetB(double b) {m_b = b;};
  double Inversa(double x) const {return (1-sqrt(1-x));};

	virtual double Eval(double x) const {return m_a*x + m_b;};

	virtual void Print() const;
	
private:
	double m_a, m_b;

};

#endif
