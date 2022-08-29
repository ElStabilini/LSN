#ifndef __RANDOM__H
#define __RANDOM__H

#ifndef __FUNCTION__H
#define __FUNCTION__H

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>

#include <vector>
#include <string>

#include "random.h"
#include "function.cpp"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;

   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

  
//DEFINIZIONE VARIABILI
  double d = 1;
  double L = (d - 0.223);
  double M = 100000; //numero di coordinate che devo estrarre (una cartesiana e una angolare)
	double Nth=500;
  int Lbox=int(M/Nth); //numero di elementi per blocco
  
  vector<double> x; //numero di lanci (ovvero numero di blocchi)
  vector<double> pi_ave; //vettore di N valori - contiene le medie di ciascun blocco
  vector<double> pi_ave2; //vettore con gli stessi elementi ma al quadrato
  vector<double> pi_sum; //questo vettore deve contenere le stime PROGRESSIVE di pi
  vector<double> pi_sum2; //contiene gli stessi elementi ma al quadrato
  vector<double> err_prog; //errore sui risultati progressivi
   
  vector<double> x1;
  vector<double> y2;
  vector<double> x2;
     
//ESTRAZIONE NUMERI CALCOLO PI 
   
  for(int i=0; i<M; i++) {
		x1.push_back(rnd.Rannyu()); //riempio r con M numeri casuali (estraggo una coordinata x - sarà il centro della circonferenza-)
		y2.push_back(rnd.Rannyu(-L,L));
	 }
	 
	for(int i=0; i<M; i++) {
		if(abs(y2[i]) == L) //se la coordinata y è (+-)L allora la sbarretta è verticale
			 x2.push_back(x1[i]);
		else {
		//direzione: prendo casualmente una direzione tra destra e sinistra (estraggo g, se g
			 double g = rnd.Rannyu(-1,1);
			 if (g>0) {
				 //prendo la coordinata destra quindi considero la cooridinata positiva
				 x2.push_back(sqrt(pow(L,2)-pow(y2[i],2))); 
				 }
			 else {
				 //prendo la coordinata sinistra
				 x2.push_back(-sqrt( pow(L,2) - pow(y2[i],2) ));
				 }
		 }
	}
	
//CALCOLO PI & STATISTICA a BLOCCHI
	BlochStatPi(Nth, L, Lbox, x1, x2, y2, pi_ave, pi_ave2, pi_sum, pi_sum2, err_prog);
	Print(pi_sum, err_prog, Nth, "Pi.dat");
	
  rnd.SaveSeed();
    
  return 0;
}

#endif
#endif
